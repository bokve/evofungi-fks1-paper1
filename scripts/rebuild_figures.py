"""Rebuild manuscript plots from fixed numeric records; no model fitting."""
from pathlib import Path
import argparse,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from paper1.figures import run
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',type=Path,default=ROOT/'data');p.add_argument('--out',type=Path,default=ROOT/'rebuilt_figures');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
inputs={str(f.relative_to(a.data_root/'results')):f for f in (a.data_root/'results').rglob('*') if f.is_file()}
run(inputs,a.out,{})
print(a.out)
