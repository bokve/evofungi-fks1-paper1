"""Single local launcher for already-authenticated, explicitly selected Colab.

Never provisions/stops a VM or collects OAuth codes. Colab CLI 0.7.4 can
return zero after a kernel exception, so a unique remote success marker is
mandatory. This is readiness infrastructure, not the scientific pipeline.
"""
import argparse
import hashlib
import json
import re
import subprocess
import tempfile
import uuid
from pathlib import Path

STAGES = {'smoke': 'remote_smoke.py', 'probe': 'remote_access_probe.py', 'intake': 'remote_intake.py', 'source_probe': 'remote_source_probe.py', 'prepare_pipeline': 'remote_prepare_pipeline.py', 'provenance_audit': 'remote_provenance_audit.py', 'provenance_followup': 'remote_provenance_followup.py', 'esm_cpu_check': 'remote_esm_cpu_check.py', 'expected_error': 'remote_expected_error.py', 'close_provenance': 'remote_close_provenance.py', 'archive_readiness': 'remote_archive_readiness.py', 'verify_delivery': 'remote_verify_delivery.py'}


def redact(text):
    """Official CLI network errors can echo short-lived proxy query tokens."""
    text=re.sub(r'([?&](?:colab-runtime-proxy-token|access_token|refresh_token|id_token|code|state)=)[^&\s\"\'<>]+',r'\1[REDACTED]',text,flags=re.I)
    return re.sub(r'(Bearer\s+)[A-Za-z0-9._~+/-]+',r'\1[REDACTED]',text,flags=re.I)


def confirmed(stdout, marker, code_hash):
    records = [json.loads(line[len(marker):]) for line in stdout.splitlines() if line.startswith(marker)]
    return len(records) == 1 and records[0] == {'success': True, 'code_sha256': code_hash}


def build_wrapper(payload, name, marker, code_hash):
    # Each CLI call shares a Jupyter kernel. A fresh process prevents stale
    # sys.modules/package paths from silently overriding the transmitted code.
    return '\n'.join([
        'import json, subprocess, sys, tempfile',
        'from pathlib import Path',
        "with tempfile.TemporaryDirectory(prefix='evofungi-stage-', dir='/tmp') as _task_dir:",
        f'    _task_path = Path(_task_dir) / {name!r}',
        f'    _task_path.write_bytes(bytes.fromhex({payload.hex()!r}))',
        '    _task_result = subprocess.run([sys.executable, "-u", str(_task_path)], capture_output=True, text=True)',
        '    if _task_result.stdout: print(_task_result.stdout, end="", flush=True)',
        '    if _task_result.stderr: print(_task_result.stderr, end="", file=sys.stderr, flush=True)',
        '    if _task_result.returncode:',
        '        raise RuntimeError(f"Isolated remote stage failed with child exit code {_task_result.returncode}")',
        f'print({marker!r} + json.dumps({{"success": True, "code_sha256": {code_hash!r}}}))',
        '',
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', required=True)
    parser.add_argument('--stage', choices=[*STAGES, 'all'], default='all')
    parser.add_argument('--script', type=Path, help='Explicit small remote entry point; runs only this script')
    parser.add_argument('--timeout', type=int, default=900, help='Remote transport timeout in seconds')
    parser.add_argument('--env', action='append', default=[], help='Non-secret remote KEY=VALUE setting')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cli = root / '.venv-tools/bin/colab'
    if not cli.is_file():
        raise SystemExit('Official Colab CLI not installed')
    if not (Path.home() / '.config/colab-cli/token.json').is_file():
        raise SystemExit('Complete official OAuth in your local terminal first; never send codes to chat')
    log_root = root / 'progress/logs'
    log_root.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='remote-readiness-', dir=log_root))
    selected = [args.script.stem] if args.script else (('smoke', 'probe', 'intake') if args.stage == 'all' else [args.stage])
    for stage in selected:
        source = args.script.resolve(strict=True) if args.script else root / 'tools' / STAGES[stage]
        payload = source.read_bytes()
        code_hash = hashlib.sha256(payload).hexdigest()
        marker = f'EVOFUNGI_REMOTE_RESULT_{uuid.uuid4().hex} '
        wrapper = build_wrapper(payload, source.name, marker, code_hash)
        with tempfile.TemporaryDirectory(prefix='evofungi-transport-') as d:
            path = Path(d) / 'stage.py'
            path.write_text(wrapper)
            command = [str(cli), '--auth', 'oauth2', 'exec', '-s', args.session, '-f', str(path), '--timeout', str(args.timeout)]
            for item in args.env:
                command.extend(['--env', item])
            print(f'Executing remote readiness {stage}; logs: {run}', flush=True)
            result = subprocess.run(command, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=args.timeout + 60)
        safe_stdout,safe_stderr=redact(result.stdout),redact(result.stderr)
        (run / f'{stage}.stdout.txt').write_text(safe_stdout)
        (run / f'{stage}.stderr.txt').write_text(safe_stderr)
        success = result.returncode == 0 and confirmed(result.stdout, marker, code_hash)
        record = {'stage': stage, 'session': args.session, 'code_sha256': code_hash, 'CLI_returncode': result.returncode, 'remote_success_marker': success, 'scientific_result': None if args.script else False, 'scope': 'transport_validation_only; scientific status belongs to stage outputs'}
        (run / f'{stage}.json').write_text(json.dumps(record, indent=2))
        print(safe_stdout, end='')
        print(safe_stderr, end='')
        print(json.dumps(record), flush=True)
        if not success:
            raise SystemExit(f'Remote {stage} failed; inspect preserved logs. No downstream stage executed.')


if __name__ == '__main__':
    main()
