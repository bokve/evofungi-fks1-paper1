"""Read-only Drive inventory; outputs stay in ephemeral Colab /content.

No notebook execution, model loading, source rewriting, or inferred schema mapping.
Large files are inventoried but explicitly left unhashed for targeted follow-up.
"""
import csv
import hashlib
import json
import os
import platform
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path('/content/drive/MyDrive/EvoFungi_AMR')
HASH_LIMIT = 16 * 1024 * 1024
TEXT_SUFFIXES = {'.csv', '.tsv', '.json', '.py', '.r', '.fasta', '.fa', '.fna', '.faa', '.md', '.txt', '.yaml', '.yml'}


def inspect_file(path, root=ROOT):
    before = path.stat()
    record = {'path': str(path.relative_to(root)), 'bytes': before.st_size, 'mtime_ns': before.st_mtime_ns, 'sha256': None, 'status': 'INVENTORIED_ONLY'}
    if path.suffix.lower() not in TEXT_SUFFIXES or before.st_size > HASH_LIMIT:
        record['reason'] = 'unselected_extension_or_over_16MiB'
        return record
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    record.update(sha256=h.hexdigest(), status='HASHED')
    try:
        if path.suffix.lower() in {'.csv', '.tsv'}:
            with path.open(encoding='utf-8-sig', newline='') as f:
                reader = csv.reader(f, delimiter='\t' if path.suffix.lower() == '.tsv' else ',')
                columns = next(reader)
                n = 0
                malformed = []
                for line, row in enumerate(reader, 2):
                    n += 1
                    if len(row) != len(columns) and len(malformed) < 10:
                        malformed.append(line)
            record.update(columns=columns, rows=n, duplicate_columns=len(columns) != len(set(columns)), malformed_row_examples=malformed)
            if record['duplicate_columns'] or malformed:
                record.update(status='SCHEMA_ERROR', error='Duplicate columns or inconsistent CSV row widths')
        elif path.suffix.lower() == '.json':
            with path.open(encoding='utf-8-sig') as f:
                obj = json.load(f)
            record['json_type'] = type(obj).__name__
            if isinstance(obj, dict):
                record['json_top_keys'] = sorted(obj)
    except (UnicodeError, csv.Error, json.JSONDecodeError, StopIteration) as e:
        record.update(status='SCHEMA_ERROR', error=f'{type(e).__name__}: {e}')
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        record.update(status='CHANGED_DURING_READ', sha256=None)
    return record


def main():
    if not ROOT.is_dir():
        raise RuntimeError('Remote project root missing; resolve account/mount before intake')
    output = Path(tempfile.mkdtemp(prefix='evofungi-intake-', dir='/content'))
    records, errors = [], []

    def walk_error(error):
        errors.append({'path': error.filename, 'error': str(error)})

    for parent, dirs, names in os.walk(ROOT, followlinks=False, onerror=walk_error):
        dirs.sort()
        for name in sorted(names):
            path = Path(parent) / name
            try:
                if path.is_symlink():
                    records.append({'path': str(path.relative_to(ROOT)), 'status': 'SYMLINK_SKIPPED'})
                else:
                    records.append(inspect_file(path))
            except OSError as e:
                errors.append({'path': str(path.relative_to(ROOT)), 'error': str(e)})
    inventory = output / 'inventory.json'
    inventory.write_text(json.dumps(records, indent=2))
    summary = {'purpose': 'read_only_remote_intake_not_scientific_validation', 'root': str(ROOT), 'python': platform.python_version(), 'inventory': str(inventory), 'inventory_sha256': hashlib.sha256(inventory.read_bytes()).hexdigest(), 'files': len(records), 'statuses': dict(Counter(r['status'] for r in records)), 'errors': errors, 'hash_limit_bytes': HASH_LIMIT, 'large_files_remain_unverified': True}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if errors or any(r['status'] in {'SCHEMA_ERROR', 'CHANGED_DURING_READ'} for r in records):
        raise RuntimeError('Intake has explicit errors; inspect inventory before selecting scientific inputs')


if __name__ == '__main__':
    main()
