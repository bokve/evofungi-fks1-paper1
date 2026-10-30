"""Pin acquired weights and repeat one real score before leaving the GPU."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .checkpoint_gate import verify_acquired_bytes
from .esm_contract import read_protein
from .esm_execution import load_checkpoint
from .esm_scoring import score_masked_position
from .stages import atomic_json

def run(inputs,attempt,config):
    source=inputs['esm-score/CHECKPOINT_ACQUISITION_MANIFEST.json']
    records=json.loads(source.read_text())
    gate=verify_acquired_bytes(source)
    model,alphabet,info=load_checkpoint(records[0]['checkpoint'])
    if info['sha256']!=records[0]['sha256']:raise ValueError('Cached model byte identity changed')
    sequence=read_protein(inputs['protein'],config['protein_record_id'])
    members=pd.read_csv(inputs['esm-score/ESM_MEMBER_SCORES.csv'])
    position=int(members.aa_position.min())
    values,_=score_masked_position(model,alphabet,sequence,position)
    selected=members[(members.checkpoint==info['checkpoint'])&(members.aa_position==position)]
    np.testing.assert_allclose([values[r.mut_aa] for r in selected.itertuples()],selected.score,rtol=1e-5,atol=1e-5)
    atomic_json(attempt/'ESM_CHECKPOINT_LOCK.json',records)
    atomic_json(attempt/'ESM_REAL_REPRODUCTION.json',{'status':'EXECUTED_VALIDATED','checkpoint_byte_gate':gate,
       'checkpoint':info['checkpoint'],'checkpoint_sha256':info['sha256'],'position':position,'actual_scores_repeated':len(selected),
       'maximum_absolute_difference':float(max(abs(values[r.mut_aa]-r.score) for r in selected.itertuples())),
       'atol':1e-5,'rtol':1e-5,'scope':'one actual cached-checkpoint position repeat; no phenotype association; no full inference rerun'})
    return ['ESM_CHECKPOINT_LOCK.json','ESM_REAL_REPRODUCTION.json']
