"""Bounded Paper 1 DAG, extending the existing hash/lock stage engine."""
import json
import os
from pathlib import Path
from .stages import atomic_json, digest, run_stage, valid_completion

STAGE_ROLES = {
 'panels':['aa','phenotype','evo2_aa','gemme','protein','classified_HS1','classified_HS2','classified_HS3','author_validation_DMS_missing_estimates','author_rescue_missing_mutants'],
 'hs2-core':['ascertainment','multiplicity','author_config',*[f'frequency_{i}' for i in range(6)]],
 'history':[],
 'lineage':['method3_long','raw_nt_scores','dna'],
 'esm-smoke':['aa','protein','esm_jobs'],
 'esm-score':['aa','protein','esm_jobs'],
 'esm-finalize':['protein'],
 'compare':[],
 'support':['h_membership','alt490','multiplicity','protein','gemme_matrix','gemme','coordinate_map','gemme_matrix_qc'],
 'figures':[],
 'writeup':[],
 'verify':['protein','visual_review'],
 'release':['protein','esm_jobs','context_archive','checkpoint_lock','visual_review','task_status_template'],
}
DEPENDENCIES = {'panels':[], 'hs2-core':[], 'history':[], 'esm-smoke':['panels'], 'esm-score':['esm-smoke']}
DEPENDENCIES['lineage']=[]
DEPENDENCIES.update({'compare':['panels','esm-score'],'support':['hs2-core','esm-score'],'figures':['panels','hs2-core','history','esm-score','compare']})
DEPENDENCIES['esm-finalize']=['esm-score']
DEPENDENCIES.update({'writeup':['panels','hs2-core','history','lineage','esm-score','compare','figures'],
 'verify':['panels','esm-score','esm-finalize','compare','figures','writeup'],
 'release':['panels','hs2-core','history','lineage','esm-smoke','esm-score','esm-finalize','compare','figures','writeup','verify']})


def execute(config_path, selected):
    config_path=Path(config_path).resolve();config=json.loads(config_path.read_text())
    root=Path(config['project_root']).resolve(strict=True)
    code_root=Path(__file__).resolve().parents[1]
    protocol=code_root/config['protocol_path']
    if digest(protocol)!=config['protocol_sha256']:
        raise ValueError('Locked new-analysis protocol hash mismatch')
    output=root/config['output_relative'];output.mkdir(parents=True,exist_ok=True)
    index_path=output/'RUN_INDEX.json'
    index=json.loads(index_path.read_text()) if index_path.exists() else {'run_id':'PAPER1_EXECUTION_2026_10_18','stages':{},'status':'RUNNING'}
    choices={'all':list(STAGE_ROLES),'cpu':['panels','hs2-core','history'],'pre-esm':['panels','hs2-core','history','esm-smoke'],
             'post-esm':['lineage','compare','support','figures','writeup'],'finalize':['verify','release']}
    stages=choices.get(selected,[selected])
    def run_one(name):
        if name not in STAGE_ROLES:raise ValueError('Stage not implemented: '+name)
        if name.startswith('esm-') and os.environ.get('EVOFUNGI_REUSE_ESM')=='1':
            from .esm_reuse import verify as verify_saved
            result=verify_saved(name,index['stages'].get(name),config,root,code_root,protocol)
            index['stages'][name]=result;atomic_json(index_path,index)
            from .task_status import update
            update(root,output,index,protocol,config)
            print(json.dumps({'stage':name,**result}),flush=True)
            return
        roles=list(STAGE_ROLES[name])
        if name=='history':
            roles=list(config['integration_roles'])+['NEW2_specification','h3_paired','method1_summary','method1_components','method2_summary','method3_summary','method3_reliability','method3_closure','structure_descriptors']
            roles+=sorted(x for x in config['inputs'] if x.startswith('history_'))
            roles+=sorted(x for x in config['inputs'] if x.startswith('context_') and x!='context_archive')
        if name=='esm-score' and 'checkpoint_lock' in config['inputs']:roles.append('checkpoint_lock')
        inputs={'protocol':protocol}
        for role in roles:
            spec=config['inputs'][role]
            base=code_root if 'context_relative_path' in spec else root
            p=(base/spec.get('context_relative_path',spec.get('relative_path'))).resolve(strict=True)
            if not p.is_relative_to(base) or digest(p)!=spec['sha256']:
                raise ValueError('Pinned scientific input mismatch: '+role)
            inputs[role]=p
        for dependency in DEPENDENCIES[name]:
            dep=index['stages'].get(dependency)
            if not dep or not dep.get('directory') or not valid_completion(dep['directory'],dep['fingerprint']):
                raise ValueError('Verified dependency required: '+dependency)
            directory=Path(dep['directory']);inputs[dependency+'_completion']=directory/'completion.json'
            for rel in json.loads((directory/'completion.json').read_text())['outputs']:
                inputs[dependency+'/'+rel]=directory/rel
        if name in ['writeup','verify','release']:
            disposition=output/'SUPPORT_TERMINAL_DISPOSITION.json'
            if not disposition.is_file():raise ValueError('Separate categorical-support terminal disposition required')
            inputs['support_disposition']=disposition
            support_record=index['stages'].get('support',{})
            if support_record.get('directory') and valid_completion(support_record['directory'],support_record.get('fingerprint')):
                directory=Path(support_record['directory'])
                inputs['support_completion']=directory/'completion.json'
                for rel in json.loads((directory/'completion.json').read_text())['outputs']:inputs['support/'+rel]=directory/rel
        from .panels import prepare
        from .context_stages import hs2_core,history
        from .esm_execution import smoke,score
        from .analysis import run as compare
        from .support_analysis import run as support
        from .figures import run as figures
        from .writeup import run as writeup
        from .verification import run as verify
        from .release import run as release
        from .esm_finalize import run as esm_finalize
        from .lineage import run as lineage
        actions={'panels':prepare,'hs2-core':hs2_core,'history':history,'lineage':lineage,'esm-smoke':smoke,'esm-score':score,'esm-finalize':esm_finalize,'compare':compare,'support':support,'figures':figures,'writeup':writeup,'verify':verify,'release':release}
        module_sets={'panels':['panels.py'],'hs2-core':['context_stages.py'],'history':['context_stages.py'],
                     'lineage':['lineage.py'],
                     'esm-smoke':['esm_execution.py','esm_scoring.py','esm_contract.py'],
                     'esm-score':['esm_execution.py','esm_scoring.py','esm_contract.py'],
                     'esm-finalize':['esm_finalize.py','checkpoint_gate.py','esm_scoring.py','esm_execution.py'],
                     'compare':['analysis.py','metrics.py','comparison.py','panels.py'],
                     'support':['support_analysis.py','metrics.py'], 'figures':['figures.py'],
                     'writeup':['writeup.py'],'verify':['verification.py','metrics.py','comparison.py','esm_execution.py','esm_scoring.py'],
                     'release':['release.py']}
        modules=module_sets[name]
        if name=='esm-score' and 'checkpoint_lock' in inputs:
            from .checkpoint_gate import verify_acquired_bytes
            verify_acquired_bytes(inputs['checkpoint_lock']);modules=[*modules,'checkpoint_gate.py']
        code={p:Path(__file__).parent/p for p in modules}
        if name=='verify':code.update({'tests/'+p.name:p for p in (code_root/'tests').glob('*.py')})
        if name=='release':
            code.update({str(p.relative_to(code_root)):p for folder in ['paper1','tools','tests','configs'] for p in (code_root/folder).glob('*') if p.is_file() and p.suffix in ['.py','.json','.md']})
            # Exclude release's own previous fingerprint and per-invocation
            # REUSED flags: neither changes the generating scientific artifacts.
            # Live invocation/task status is delivered separately after snapshot.
            snapshot={'run_id':index['run_id'],'scope':'generating artifacts before release; live invocation status is separate','stages':{}}
            for stage,record in index['stages'].items():
                if stage=='release':continue
                if record.get('directory'):
                    snapshot['stages'][stage]={key:record[key] for key in ['directory','fingerprint']}
                    snapshot['stages'][stage]['completion_sha256']=digest(Path(record['directory'])/'completion.json')
                else:snapshot['stages'][stage]=record
            snapshot_path=output/'GENERATING_ARTIFACT_INDEX.json';atomic_json(snapshot_path,snapshot)
            task_snapshot_path=output/'TASK_ARTIFACT_INDEX.json'
            tasks=json.loads((output/'TASK_STATUS.json').read_text())
            task_snapshot={'scope':'Artifact disposition before release; T21/T22 and actual resume flags are resolved in post-snapshot live task status','tasks':[]}
            for task in tasks['tasks']:
                if task['id'] in ['T21','T22']:continue
                keys=['id','title','evidence_status','stage','fingerprint','completion_manifest_path','completion_manifest_sha256','block_reason','scope','actual_output_paths']
                entry={k:task[k] for k in keys if k in task}
                entry['artifact_status']='EXECUTED_VALIDATED' if task['execution_status']=='REUSED_VALIDATED' else task['execution_status']
                task_snapshot['tasks'].append(entry)
            atomic_json(task_snapshot_path,task_snapshot)
            inputs['run_index_before_release']=snapshot_path
            inputs['task_status_before_release']=task_snapshot_path
        def validate(attempt, outputs):
            for rel in outputs:
                p=attempt/rel
                empty_stream=rel.endswith(('.stdout.txt','.stderr.txt'))
                if not p.is_file() or (p.stat().st_size==0 and not empty_stream):raise ValueError('Missing/empty output: '+rel)
                if p.suffix=='.json':json.loads(p.read_text())
        try:
            stage_params={k:v for k,v in config.items() if k not in ['inputs','project_root','output_relative','protocol_path']}
            stage_params['stage']=name
            result=run_stage(name,output,inputs,code,stage_params,
                             lambda attempt:actions[name](inputs,attempt,config),validate,
                             packages=['numpy','pandas','torch','fair-esm'] if name.startswith('esm-') else ['numpy','pandas','scipy','matplotlib'])
            index['stages'][name]=result;atomic_json(index_path,index)
            if name=='support':
                gate=json.loads((Path(result['directory'])/'HS2_SUPPORT_SCORE_GATE.json').read_text())
                atomic_json(output/'SUPPORT_TERMINAL_DISPOSITION.json',{'execution_status':gate['status'],'gate_path':str(Path(result['directory'])/'HS2_SUPPORT_SCORE_GATE.json'),'scope':'categorical H208, independent from continuous HS2'})
            if name=='release':
                index['status']='ANALYSES_COMPLETE_WITH_BLOCKED_BRANCHES';atomic_json(index_path,index)
            from .task_status import update
            update(root,output,index,protocol,config)
            print(json.dumps({'stage':name,**result}),flush=True)
        except BaseException:
            index['status']='PARTIAL';atomic_json(index_path,index);raise
    failed=[]
    for name in stages:
        try:run_one(name)
        except BaseException as error:
            import traceback
            failure={'status':'BLOCKED_INPUT' if isinstance(error,(ValueError,FileNotFoundError,KeyError)) else 'FAILED_TECHNICAL','error_type':type(error).__name__,'error':str(error),'traceback':traceback.format_exc()}
            index['stages'][name]=failure;index['status']='PARTIAL';atomic_json(index_path,index)
            atomic_json(output/(name+'_latest_failure.json'),failure)
            if name=='support':atomic_json(output/'SUPPORT_TERMINAL_DISPOSITION.json',{'execution_status':failure['status'],'block_reason':failure['error'],'scope':'affected categorical feature/identity branch stopped'})
            from .task_status import update
            update(root,output,index,protocol,config)
            failed.append(name);print(json.dumps({'stage':name,**failure}),flush=True)
            if selected not in choices:raise
    print(json.dumps({'run_index':str(index_path),'status':index['status'],'failed_stages':failed}),flush=True)
    if any(name!='support' for name in failed):raise RuntimeError('Required/dependent stages did not complete: '+','.join(failed))
