"""Real five-checkpoint ESM inference, exclusively on a remote CUDA runtime."""
import argparse
import gc
import hashlib
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
from .esm_contract import CHECKPOINTS,CANONICAL_AA,read_protein,targets,window
from .esm_scoring import masked_tokens,score_masked_position
from .stages import atomic_json,digest,run_stage,valid_completion


def environment():
    import torch
    if not torch.cuda.is_available():raise RuntimeError('Actual CUDA GPU required; never substitute CPU inference')
    if not Path('/content').is_dir():raise RuntimeError('ESM weights must stay in Colab')
    torch.manual_seed(20261018);np.random.seed(20261018)
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.backends.cudnn.deterministic=True
    torch.set_float32_matmul_precision('highest')
    torch.hub.set_dir('/content/evofungi-model-cache')
    # Official ESM-1v checkpoints contain argparse.Namespace metadata. Keep
    # weights-only loading with this narrow allowlist on current PyTorch.
    torch.serialization.add_safe_globals([argparse.Namespace])
    return torch


def load_checkpoint(name):
    import esm
    torch=environment()
    if name not in CHECKPOINTS:raise ValueError('Unfrozen checkpoint name')
    model,alphabet=esm.pretrained.load_model_and_alphabet(name)
    model=model.eval().float().cuda()
    path=Path(torch.hub.get_dir())/'checkpoints'/(name+'.pt')
    if not path.is_file():raise RuntimeError('Downloaded checkpoint locator missing')
    info={'checkpoint':name,'acquisition_url':'https://dl.fbaipublicfiles.com/fair-esm/models/'+name+'.pt',
          'remote_cache_path':str(path),'bytes':path.stat().st_size,'sha256':digest(path),
          'dtype':str(next(model.parameters()).dtype),'device':str(next(model.parameters()).device),
          'GPU':torch.cuda.get_device_name(0),'token_capacity':getattr(model.args,'max_positions',None),
          'eval':not model.training,'scoring':'masked_marginal_ALT_minus_WT','seed':20261018}
    if info['token_capacity'] is None or info['token_capacity']<1024:
        raise ValueError('Actual model capacity incompatible; technical amendment required before inference')
    return model,alphabet,info


def verify_jobs(inputs,config):
    sequence=read_protein(inputs['protein'],config['protein_record_id'])
    aa=pd.read_csv(inputs['aa'])
    jobs=targets(sequence,aa.to_dict('records'))
    frozen=json.loads(inputs['esm_jobs'].read_text())
    if jobs!=frozen:raise ValueError('17 October frozen ESM windows/positions changed')
    if len(jobs)!=140 or len(sequence)!=1876:raise ValueError('Frozen ESM design mismatch')
    return sequence,aa,jobs


def smoke(inputs,attempt,config):
    torch=environment();sequence,aa,jobs=verify_jobs(inputs,config)
    model,alphabet,info=load_checkpoint(CHECKPOINTS[0]);torch.cuda.reset_peak_memory_stats()
    first=jobs[0]['aa_position']
    a,_=score_masked_position(model,alphabet,sequence,first)
    b,_=score_masked_position(model,alphabet,sequence,first)
    np.testing.assert_allclose(list(a.values()),list(b.values()),rtol=1e-5,atol=1e-5)
    short=sequence[:64];pos=32
    actual,_=score_masked_position(model,alphabet,short,pos)
    tokens,idx,wt,_=masked_tokens(alphabet,short,pos)
    with torch.inference_mode():
        # Official predict.py masked-marginals + label_row formula, evaluated
        # at precisely the same single mask/BOS index and actual weights.
        official=torch.log_softmax(model(tokens.cuda())['logits'],dim=-1)
    expected={alt:float((official[0,idx,alphabet.get_idx(alt)]-official[0,idx,alphabet.get_idx(wt)]).cpu()) for alt in sorted(CANONICAL_AA)}
    np.testing.assert_allclose(list(actual.values()),list(expected.values()),rtol=1e-5,atol=1e-5)
    boundaries=[]
    for p in sorted({1,len(sequence)//2,len(sequence),*[j['aa_position'] for j in jobs]}):
        tok,idx,wt,w=masked_tokens(alphabet,sequence,p)
        if sequence[p-1]!=w['sequence'][idx-1] or tok.shape[1]>info['token_capacity']:
            raise ValueError('Actual alphabet boundary/index mismatch')
        boundaries.append({'position':p,'tokens':int(tok.shape[1]),'target_token':int(idx),**{k:v for k,v in w.items() if k!='sequence'}})
    report={'status':'EXECUTED_VALIDATED','actual_weights_loaded':True,'checkpoint':info,
            'repeat_max_absolute_difference':max(abs(a[k]-b[k]) for k in a),
            'official_formula_max_absolute_difference':max(abs(actual[k]-expected[k]) for k in actual),
            'tolerance_absolute':1e-5,'tolerance_relative':1e-5,'finite_canonical_AA':True,'WT_self_zero':True,
            'boundary_token_checks':boundaries,'max_GPU_memory_bytes':torch.cuda.max_memory_allocated(),
            'phenotype_association_computed':False,'determinism_caveat':'Same-runtime numerical repeat check, not cross-hardware bitwise reproducibility.'}
    atomic_json(attempt/'ESM_GPU_SMOKE.json',report)
    atomic_json(attempt/'CHECKPOINT_ACQUISITION_MANIFEST.json',[info])
    del model;gc.collect();torch.cuda.empty_cache()
    return ['ESM_GPU_SMOKE.json','CHECKPOINT_ACQUISITION_MANIFEST.json']


def score(inputs,attempt,config):
    torch=environment();sequence,aa,jobs=verify_jobs(inputs,config)
    smoke_report=json.loads(inputs['esm-smoke/ESM_GPU_SMOKE.json'].read_text())
    if smoke_report['status']!='EXECUTED_VALIDATED' or not smoke_report['actual_weights_loaded']:
        raise ValueError('Real checkpoint smoke required')
    output=Path(config['project_root'])/config['output_relative']/'esm_position_jobs'
    tables=[];acquisition=[];job_records=[]
    for name in CHECKPOINTS:
        model,alphabet,info=load_checkpoint(name);acquisition.append(info)
        if name==CHECKPOINTS[0] and info['sha256']!=smoke_report['checkpoint']['sha256']:
            raise ValueError('Checkpoint bytes differ from real smoke')
        for job in [j for j in jobs if j['checkpoint']==name]:
            position=job['aa_position']
            def action(directory):
                start=time.monotonic();values,w=score_masked_position(model,alphabet,sequence,position)
                atomic_json(directory/'scores.json',{'checkpoint':name,'aa_position':position,'wt_aa':job['wt_aa'],
                     'scores':values,'window':w,'checkpoint_sha256':info['sha256'],'wall_seconds':time.monotonic()-start})
                return ['scores.json']
            def validate(directory,outputs):
                r=json.loads((directory/'scores.json').read_text())
                if set(r['scores'])!=CANONICAL_AA or r['scores'][job['wt_aa']]!=0 or not all(np.isfinite(list(r['scores'].values()))):
                    raise ValueError('Incomplete or nonfinite real ESM scores')
            code={n:Path(__file__).parent/n for n in ['esm_execution.py','esm_scoring.py','esm_contract.py']}
            result=run_stage('member_'+name+'_position_'+str(position),output,
                   {'protein':inputs['protein'],'protocol':inputs['protocol'],'frozen_jobs':inputs['esm_jobs']},code,
                   {'checkpoint':info,'job':job,'precision':'float32','seed':20261018},action,validate,packages=['torch','fair-esm','numpy'])
            r=json.loads((Path(result['directory'])/'scores.json').read_text());job_records.append(result)
            for alt,val in r['scores'].items():
                if alt==job['wt_aa']:continue
                tables.append({'aa_variant_id':f"FKS1:{job['wt_aa']}{position}{alt}",'aa_position':position,'wt_aa':job['wt_aa'],'mut_aa':alt,
                     'checkpoint':name,'score':val,'checkpoint_sha256':info['sha256'],'window_sha256':job['window_sha256'],
                     'start_0based':job['start_0based'],'end_0based_exclusive':job['end_0based_exclusive'],
                     'measured_RAW435':f"FKS1:{job['wt_aa']}{position}{alt}" in set(aa.aa_variant_id)})
            print(json.dumps({'completed_ESM_position':position,'checkpoint':name,'status':result['status']}),flush=True)
        path=Path(info['remote_cache_path'])
        if digest(path)!=info['sha256']:raise ValueError('Checkpoint changed during scoring')
        del model;gc.collect();torch.cuda.empty_cache()
        atomic_json(attempt/'CHECKPOINT_ACQUISITION_MANIFEST.json',acquisition)
    members=pd.DataFrame(tables)
    if len(members)!=5*532 or members.duplicated(['aa_variant_id','checkpoint']).any():raise ValueError('Five-member exact coverage failed')
    pivot=members.pivot(index='aa_variant_id',columns='checkpoint',values='score').reindex(columns=CHECKPOINTS)
    if pivot.shape!=(532,5) or not np.isfinite(pivot.to_numpy()).all():raise ValueError('No partial ensemble allowed')
    ensemble=members.drop_duplicates('aa_variant_id')[['aa_variant_id','aa_position','wt_aa','mut_aa','measured_RAW435']].set_index('aa_variant_id')
    ensemble['esm1v_ensemble']=pivot.mean(axis=1)
    if set(aa.aa_variant_id)-set(ensemble.index):raise ValueError('Missing primary AA ESM coverage')
    members.to_csv(attempt/'ESM_MEMBER_SCORES.csv',index=False)
    ensemble.reset_index().to_csv(attempt/'ESM_ENSEMBLE.csv',index=False)
    ensemble.loc[aa.aa_variant_id].reset_index().to_csv(attempt/'ESM_RAW435.csv',index=False)
    atomic_json(attempt/'ESM_COVERAGE.json',{'checkpoints':5,'position_jobs':len(job_records),'positions':28,'primary_AA':435,'theoretical_AA_no_phenotype':532,'all_finite':True,'partial_ensemble':False})
    atomic_json(attempt/'ESM_JOB_INDEX.json',job_records)
    return ['ESM_MEMBER_SCORES.csv','ESM_ENSEMBLE.csv','ESM_RAW435.csv','ESM_COVERAGE.json','ESM_JOB_INDEX.json','CHECKPOINT_ACQUISITION_MANIFEST.json']
