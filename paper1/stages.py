"""Outcome-independent, hash-verified stage storage for CPU/GPU work."""
import hashlib
import importlib.metadata
import json
import os
import platform
import traceback
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, record):
    path = Path(path)
    temporary = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    try:
        with temporary.open('x') as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def writer_lock(output_root):
    """Exclusive mkdir lock; never silently break another/stale writer lock.

    Tested for processes sharing one filesystem view. Multi-runtime Drive
    FUSE consistency is not guaranteed; use one designated Colab writer.
    """
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    lock = root / '.paper1-writer.lock'
    owner = uuid.uuid4().hex
    try:
        lock.mkdir()
    except FileExistsError:
        raise RuntimeError(f'Writer lock exists: {lock}. Resolve the owning process; no automatic stale-lock deletion.') from None
    try:
        atomic_json(lock / 'owner.json', {'owner': owner, 'pid': os.getpid(), 'host': platform.node(), 'created_utc': datetime.now(timezone.utc).isoformat()})
        yield
    finally:
        owner_file = lock / 'owner.json'
        if owner_file.is_file() and json.loads(owner_file.read_text()).get('owner') == owner:
            owner_file.unlink()
            lock.rmdir()


def identity(paths):
    records = {}
    for label, path in sorted(paths.items()):
        p = Path(path).resolve(strict=True)
        if not p.is_file():
            raise ValueError(f'Input is not a file: {p}')
        records[label] = {'path': str(p), 'bytes': p.stat().st_size, 'sha256': digest(p)}
    return records


def software(packages=()):
    return {'python': platform.python_version(), 'platform': platform.platform(), 'packages': {name: importlib.metadata.version(name) for name in sorted(packages)}}


def valid_completion(directory, fingerprint):
    directory = Path(directory)
    try:
        completion = json.loads((directory / 'completion.json').read_text())
        if completion['fingerprint'] != fingerprint or completion['status'] != 'EXECUTED_VALIDATED' or not completion['outputs']:
            return False
        for relative, record in completion['outputs'].items():
            p = (directory / relative).resolve(strict=True)
            if not p.is_relative_to(directory.resolve()) or not p.is_file() or p.stat().st_size != record['bytes'] or digest(p) != record['sha256']:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_stage(name, output_root, inputs, code, params, action, validate, packages=()):
    """Run only an explicitly implemented action, then validate before completion.

    Exceptions are preserved and propagated; no automatic retry of data errors.
    `action(attempt)` writes new files and returns their relative paths.
    """
    if not name.replace('_', '').replace('-', '').isalnum():
        raise ValueError('Invalid stage name')
    code = dict(code)
    code['stage_engine'] = Path(__file__)
    before = {'inputs': identity(inputs), 'code': identity(code), 'params': params, 'software': software(packages)}
    # Absolute deployment paths are provenance, not content identity. Otherwise
    # an identical bundle extracted in a new /content directory defeats resume.
    content_identity = dict(before)
    for field in ['inputs', 'code']:
        content_identity[field] = {key: {k:v for k,v in value.items() if k != 'path'}
                                   for key,value in before[field].items()}
    fingerprint = hashlib.sha256(json.dumps(content_identity, sort_keys=True, allow_nan=False).encode()).hexdigest()
    with writer_lock(output_root):
        base = Path(output_root) / name / fingerprint
        base.mkdir(parents=True, exist_ok=True)
        for attempt in sorted(base.iterdir()):
            if attempt.is_dir() and valid_completion(attempt, fingerprint):
                return {'status': 'REUSED_VALIDATED', 'directory': str(attempt), 'fingerprint': fingerprint}
        attempt = base / ('attempt-' + uuid.uuid4().hex)
        attempt.mkdir()
        atomic_json(attempt / 'started.json', dict(before, fingerprint=fingerprint, stage=name, status='RUNNING'))
        try:
            outputs = list(action(attempt))
            if not outputs or len(outputs) != len(set(outputs)):
                raise ValueError('Stage must declare nonempty, unique output files')
            for relative in outputs:
                p = attempt / relative
                if p.is_symlink() or not p.resolve(strict=True).is_relative_to(attempt.resolve()):
                    raise ValueError('Output must be a regular file inside this attempt')
            validate(attempt, outputs)
            after = {'inputs': identity(inputs), 'code': identity(code), 'params': params, 'software': software(packages)}
            if before != after:
                raise RuntimeError('Input/code/software changed during execution; no completion written')
            hashes = {relative: {'sha256': digest(attempt / relative), 'bytes': (attempt / relative).stat().st_size} for relative in outputs}
            atomic_json(attempt / 'completion.json', dict(before, fingerprint=fingerprint, stage=name, status='EXECUTED_VALIDATED', completed_utc=datetime.now(timezone.utc).isoformat(), outputs=hashes))
            return {'status': 'EXECUTED_VALIDATED', 'directory': str(attempt), 'fingerprint': fingerprint}
        except BaseException as error:
            atomic_json(attempt / 'failure.json', {'status': 'FAILED_TECHNICAL' if not isinstance(error, (ValueError, FileNotFoundError)) else 'BLOCKED_INPUT', 'error_type': type(error).__name__, 'error': str(error), 'traceback': traceback.format_exc()})
            raise
