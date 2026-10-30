"""Acquired ESM checkpoint-byte pin for subsequent clean/resumed execution."""
import json
from pathlib import Path
from .esm_contract import CHECKPOINTS
from .stages import digest

def verify_acquired_bytes(lock_path):
    import torch
    records=json.loads(Path(lock_path).read_text())
    if {r['checkpoint'] for r in records}!=set(CHECKPOINTS) or len(records)!=5:
        raise ValueError('All five acquired checkpoint identities must be pinned')
    torch.hub.set_dir('/content/evofungi-model-cache')
    directory=Path(torch.hub.get_dir())/'checkpoints';directory.mkdir(parents=True,exist_ok=True)
    for row in records:
        name=row['checkpoint'];p=directory/(name+'.pt')
        expected_url='https://dl.fbaipublicfiles.com/fair-esm/models/'+name+'.pt'
        if row['acquisition_url']!=expected_url:raise ValueError('Noncanonical ESM acquisition locator')
        if not p.is_file():
            temporary=p.with_suffix('.pt.partial')
            torch.hub.download_url_to_file(expected_url,str(temporary),progress=False)
            if digest(temporary)!=row['sha256'] or temporary.stat().st_size!=row['bytes']:
                raise ValueError('Acquired checkpoint bytes changed; no new inference or substitution allowed')
            temporary.replace(p)
        elif digest(p)!=row['sha256'] or p.stat().st_size!=row['bytes']:
            raise ValueError('Cached checkpoint bytes differ from acquired lock')
    return {'status':'EXECUTED_VALIDATED','checkpoint_count':5,'lock_sha256':digest(lock_path)}
