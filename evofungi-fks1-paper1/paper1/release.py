"""Small versioned scientific release with actual Drive archive readback."""
import hashlib
import json
import os
import shutil
import zipfile
from pathlib import Path
from .stages import atomic_json,digest


def run(inputs,attempt,config):
    code_root=Path(__file__).resolve().parents[1]
    files={}
    for role,path in inputs.items():
        if role=='protocol':files['protocol/MACHINE_READABLE_PROTOCOL.json']=path
        elif '/' in role:
            files['results/'+role]=path
        elif role in ['context_archive','visual_review','checkpoint_lock','task_status_template','run_index_before_release','task_status_before_release','support_disposition']:
            files['context/'+path.name]=path
    for name in ['paper1','tests','tools','configs','progress/execution_2026-10-18']:
        for p in (code_root/name).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':files['code/'+str(p.relative_to(code_root))]=p
    for name in ['AGENTS.md','CURRENT_STATE.md']:
        p=code_root/name
        if p.is_file():files['live_context_at_bundle/'+name]=p
    # Preserve the actual launchable directory, including real test fixtures
    # and explicit historical-context fallbacks used by the transport launcher.
    transport_manifest=json.loads((code_root/'code_manifest.json').read_text())
    for relative,expected in transport_manifest.items():
        p=(code_root/relative).resolve(strict=True)
        if not p.is_relative_to(code_root) or digest(p)!=expected:raise ValueError('Release code source differs from verified bundle: '+relative)
        files['code/'+relative]=p
    files['code/code_manifest.json']=code_root/'code_manifest.json'
    # Include small verified reference/target inputs; large raw data and weights
    # are represented by exact locators+hashes, not duplicated in the release.
    for role in ['protein','esm_jobs']:
        if role in inputs:files['small_inputs/'+inputs[role].name]=inputs[role]
    locators={role:spec for role,spec in config['inputs'].items()}
    atomic_json(attempt/'INPUT_LOCATOR_MANIFEST.json',{'project_root':config['project_root'],'roles':locators,'scope':'authenticated original Drive inputs; validated context fallbacks explicit'})
    files['INPUT_LOCATOR_MANIFEST.json']=attempt/'INPUT_LOCATOR_MANIFEST.json'
    software={role:json.loads(path.read_text())['software'] for role,path in inputs.items()
              if role.endswith('_completion')}
    import importlib.metadata
    import subprocess
    r_version=subprocess.run(['Rscript','--version'],capture_output=True,text=True,check=True)
    software['release_export_tools']={'python-docx':importlib.metadata.version('python-docx'),
                                      'R':(r_version.stdout or r_version.stderr).strip()}
    atomic_json(attempt/'SOFTWARE_LOCK.json',{'generating_stage_environments':software,
        'scope':'Actual Python/platform and relevant package versions for each completed stage; GPU generating environment retained under CPU reuse. Hosted historical installed environment is not reconstructed.'})
    files['SOFTWARE_LOCK.json']=attempt/'SOFTWARE_LOCK.json'
    readme='''# EvoFungi-AMR Paper 1 computational release

Scope: bounded Fks1-only execution. Historical modules are reused records, not newly fitted results. New ESM inference and comparisons are separated by source/status. Continuous PRE→T0 remains BLOCKED_INPUT; NEW-3I and the optional separate same-Fks1 validation comparison remain DEFERRED_SCOPE. No clinical/lethality/equivalence or full-study preregistration claim.

Actual resume/reproduce command from the code directory (official CLI authentication and the selected mounted Colab CPU session required):

    python3 tools/run_paper1.py --session evofungi-paper1 --stage all --reuse-esm

This command verifies/reuses saved GPU inference and reproduces downstream work on the selected CPU. It checks the generating code, inputs, protocol, checkpoint identities and persisted output hashes; it does not pretend the CPU performed model inference. For a fresh full inference use the same launcher with --session CPU_NAME --gpu-session GPU_NAME --stage all, with official Drive mounts on both sessions. The sessions are ephemeral; create the needed authorized session and complete official browser consent if expired. Preserve unrelated runtimes and run only one scientific Drive writer. The launcher transfers small code/context; raw scientific data and persisted outputs remain on Drive, and fresh model acquisition occurs only on Colab. The completed task GPU was released, so its temporary weight cache is not assumed to persist; saved scores and acquired file hashes remain available. Real smoke and all five checkpoint identities are mandatory. Acquired byte locks govern subsequent model acquisition. No historical notebook run-all.

A fresh local code checkout needs the official CLI at .venv-tools/bin/colab (tested google-colab-cli==0.7.4); create that local tooling environment and use the official browser OAuth flow. The virtual environment and credentials are deliberately excluded. A fresh remote CPU needs the recorded NumPy/Pandas/SciPy/Matplotlib/fair-esm versions; GPU inference additionally needs the generating PyTorch/CUDA-compatible environment recorded in SOFTWARE_LOCK.json. The command assumes the authenticated original Drive project and pinned input locators remain accessible.

Source-complete context is included as a separately verified ZIP under context; its packaged 17 October state is historical. The latest live state is delivered separately as post-snapshot metadata to avoid self-referential hashes. Input locators/hashes describe scientific data not duplicated here. Original H3 NPZ arrays and hosted API request bytes are not claimed freshly reconstructed. The immutable hosted Evo2 checkpoint SHA remains UNVERIFIED. Giant ESM weights, caches, virtual environments and authorization state are excluded. Unknown author administrative fields require author input.
'''
    (attempt/'README.md').write_text(readme);files['README.md']=attempt/'README.md'
    manifest={'files':{name:{'bytes':p.stat().st_size,'sha256':digest(p)} for name,p in sorted(files.items())},
       'scope':'code, acquired checkpoint identities, small input locators, actual outputs/source tables, figures/drafts/verification and source-complete context',
       'excluded':['credentials','browser/OAuth state','virtualenvs','model weights/cache','unnecessary large raw data'],
       'known_limits':['hosted Evo2 immutable checkpoint SHA unverified','historical H3 arrays and hosted requests not freshly reconstructed','continuous HS2 input blocked']}
    atomic_json(attempt/'RELEASE_MANIFEST.json',manifest)
    archive=attempt/'release.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,p in sorted(files.items()):z.write(p,name)
        z.write(attempt/'RELEASE_MANIFEST.json','RELEASE_MANIFEST.json')
    import tempfile
    import sys
    with tempfile.TemporaryDirectory(prefix='paper1-release-restore-',dir='/content') as temporary:
        restored=Path(temporary)
        with zipfile.ZipFile(archive) as z:
            for name in z.namelist():
                if name.startswith('code/'):
                    if not (restored/name).resolve().is_relative_to(restored):raise ValueError('Unsafe release code member')
                    z.extract(name,restored)
        restored_tests=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=restored/'code',env=dict(os.environ,PYTHONPATH=str(restored/'code')),capture_output=True,text=True)
        (attempt/'release_restore_tests.stdout.txt').write_text(restored_tests.stdout)
        (attempt/'release_restore_tests.stderr.txt').write_text(restored_tests.stderr)
        if restored_tests.returncode:raise RuntimeError('Tests from actual release-restored code failed; preserve archive as unsuccessful attempt')
    sha=digest(archive)
    destination=Path(config['project_root'])/'07_reports/PAPER1_EXECUTION_2026_10_18/releases'/sha
    destination.mkdir(parents=True,exist_ok=False)
    remote=destination/'release.zip';shutil.copy2(archive,remote)
    if digest(remote)!=sha:raise ValueError('Remote archive readback SHA mismatch')
    with zipfile.ZipFile(remote) as z:
        actual=json.loads(z.read('RELEASE_MANIFEST.json'))
        if actual!=manifest:raise ValueError('Remote manifest readback differs')
        for name,row in actual['files'].items():
            data=z.read(name)
            if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:raise ValueError('Remote archive member mismatch: '+name)
    shutil.copy2(attempt/'RELEASE_MANIFEST.json',destination/'RELEASE_MANIFEST.json')
    required=['esm-score/ESM_COVERAGE.json','compare/COMPARISON_SUMMARIES.csv','writeup/WRITEUP_STATUS.json','verify/REPRODUCIBILITY_REPORT.json','verify/FIGURE_INDEX_VERIFIED.json']
    if any(key not in inputs for key in required):raise ValueError('Mandatory scientific/writeup/verification outputs required')
    blocked=['HS2 continuous PRE-to-T0: comparable units/count-denominator provenance unavailable']
    support=json.loads(inputs['support_disposition'].read_text())
    if support['execution_status']!='EXECUTED_VALIDATED':blocked.append('Categorical H208 score-support: '+support.get('block_reason',support['execution_status']))
    report={'status':'RELEASE_VERIFIED','remote_directory':str(destination),'remote_archive':str(remote),'sha256':sha,'bytes':remote.stat().st_size,
         'members_verified':len(manifest['files']),'analysis_status':'ANALYSES_COMPLETE_WITH_BLOCKED_BRANCHES',
         'release_restored_code_tests_returncode':restored_tests.returncode,
         'blocked_branches':blocked,
         'manuscript_status':'MANUSCRIPT_DRAFT_COMPLETE','post_snapshot_state_delivery':'separate delivery index required; do not overwrite this archive',
         'resume_command':'python3 tools/run_paper1.py --session evofungi-paper1 --stage all --reuse-esm'}
    atomic_json(attempt/'RELEASE_DELIVERY.json',report)
    return ['INPUT_LOCATOR_MANIFEST.json','SOFTWARE_LOCK.json','README.md','RELEASE_MANIFEST.json','release.zip','RELEASE_DELIVERY.json','release_restore_tests.stdout.txt','release_restore_tests.stderr.txt']
