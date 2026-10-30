"""Verify saved GPU inference before reproducing downstream work on CPU."""
import json
from pathlib import Path
from .stages import digest,valid_completion
from .esm_contract import CHECKPOINTS

def verify(name,record,config,root,code_root,protocol):
    if not record or not record.get('directory') or not valid_completion(record['directory'],record.get('fingerprint')):
        raise ValueError('Complete saved GPU output required for --reuse-esm: '+name)
    directory=Path(record['directory']);completion=json.loads((directory/'completion.json').read_text())
    for role in ['aa','protein','esm_jobs']:
        if role not in completion['inputs']:continue
        spec=config['inputs'][role];p=Path(root)/spec['relative_path']
        if digest(p)!=spec['sha256'] or completion['inputs'][role]['sha256']!=spec['sha256']:
            raise ValueError('Saved GPU input changed: '+role)
    if digest(protocol)!=completion['inputs']['protocol']['sha256'] or completion['params']['protocol_sha256']!=config['protocol_sha256']:
        raise ValueError('Saved GPU protocol differs; no inference reuse')
    if completion['params'].get('protein_record_id')!=config['protein_record_id']:
        raise ValueError('Saved protein record selection differs')
    for label,spec in completion['code'].items():
        source=Path(__file__).parent/('stages.py' if label=='stage_engine' else label)
        if not source.is_file() or digest(source)!=spec['sha256']:
            raise ValueError('Saved inference code changed; explicit GPU rerun required: '+label)
    if completion['software']['packages'].get('fair-esm')!='2.0.0':raise ValueError('Unexpected recorded inference software')
    if name=='esm-score':
        coverage=json.loads((directory/'ESM_COVERAGE.json').read_text())
        if any(coverage[k]!=v for k,v in [('checkpoints',5),('position_jobs',140),('positions',28),('primary_AA',435),('theoretical_AA_no_phenotype',532)]):
            raise ValueError('Saved ensemble coverage incomplete')
        members=json.loads((directory/'CHECKPOINT_ACQUISITION_MANIFEST.json').read_text())
        if len(members)!=5 or {r['checkpoint'] for r in members}!=set(CHECKPOINTS):raise ValueError('Saved checkpoint identity set differs')
        if 'checkpoint_lock' in config['inputs']:
            spec=config['inputs']['checkpoint_lock'];base=Path(code_root) if 'context_relative_path' in spec else Path(root)
            p=base/spec.get('context_relative_path',spec.get('relative_path'))
            if digest(p)!=spec['sha256'] or json.loads(p.read_text())!=members:raise ValueError('Acquired checkpoint byte lock differs')
    return dict(record,status='REUSED_VALIDATED',reuse_scope='saved GPU inference bytes/code/inputs/protocol verified; generating GPU software retained, no CPU inference substitute')
