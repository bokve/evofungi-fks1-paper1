"""Check file identities without executing analysis code."""
from pathlib import Path
import argparse,json,hashlib
p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path,nargs='?',default=Path(__file__).resolve().parents[1]);a=p.parse_args();m=json.loads((a.root/'FILE_MANIFEST.json').read_text());bad=[]
for name,r in m['files'].items():
 f=a.root/name
 if not f.is_file() or f.stat().st_size!=r['bytes'] or hashlib.sha256(f.read_bytes()).hexdigest()!=r['sha256']:bad.append(name)
if bad:raise SystemExit('FAIL: '+repr(bad))
print('PASS:',len(m['files']),'file identities; this does not establish scientific validity')
