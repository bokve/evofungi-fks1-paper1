"""Run-specific task evidence/execution table, independent of planning snapshots."""
import json
from pathlib import Path
from .stages import atomic_json,digest

MAPPING={'T02':'panels','T03':'panels','T04':'panels','T05':'history','T06':'hs2-core','T08':'esm-smoke','T09':'esm-score',
 **{f'T{i:02}':'compare' for i in range(10,16)},'T16':'figures','T17':'writeup','T18':'verify','T19':'writeup','T20':'verify','T21':'release'}

def update(root,output,index,protocol,config):
    spec=config['inputs']['task_status_template'];source=Path(root)/spec['relative_path']
    if digest(source)!=spec['sha256']:raise ValueError('Historical planning matrix identity changed')
    tasks=json.loads(source.read_text());tasks['run_id']=index['run_id'];tasks['run_status']=index['status']
    for task in tasks['tasks']:
        tid=task['id'];task['evidence_status']=task['evidence_status_at_package_creation']
        stage=MAPPING.get(tid);record=index['stages'].get(stage)
        if tid=='T00':task.update(execution_status='EXECUTED_VALIDATED',actual_output_paths=[str(Path(root)/config['inputs']['context_archive']['relative_path'])],scope='340-member source-complete context readback verified; live state retained')
        elif tid=='T01':task.update(execution_status='EXECUTED_VALIDATED',actual_output_paths=[str(protocol)],protocol_sha256=config['protocol_sha256'],scope='locked before new ESM associations, historical results known')
        elif tid=='T07':
            core=index['stages'].get('hs2-core',{});path=str(Path(core['directory'])/'HS2_FEASIBILITY_GATE.json') if core.get('directory') else ''
            task.update(execution_status='BLOCKED_INPUT' if path else 'NOT_STARTED',actual_output_paths=[path] if path else [],block_reason='Continuous PRE-to-T0 comparable units/count-denominator provenance not established' if path else None)
        elif tid=='T23':task.update(execution_status='DEFERRED_SCOPE',block_reason='No separate optional same-Fks1 validation comparison; author validation source used for required replacements, not independent validation')
        elif tid=='T24':
            p=Path(output)/'SUPPORT_TERMINAL_DISPOSITION.json'
            if p.is_file():task.update(json.loads(p.read_text()));task['actual_output_paths']=[str(p)];task['disposition_sha256']=digest(p)
        elif tid=='T22':
            delivery=Path(output)/'FINAL_DELIVERY_INDEX.json'
            release=index['stages'].get('release',{})
            if delivery.is_file() and release.get('directory'):
                record=json.loads(delivery.read_text())
                release_record=json.loads((Path(release['directory'])/'RELEASE_DELIVERY.json').read_text())
                if record.get('status')=='DELIVERED_VERIFIED' and record.get('release_sha256')==release_record['sha256']:
                    for name,spec in record['files'].items():
                        if digest(Path(record['metadata_directory'])/name)!=spec['sha256']:raise ValueError('Final delivery metadata changed: '+name)
                    task.update(execution_status='EXECUTED_VALIDATED',actual_output_paths=[str(delivery)],scope='final state/report/changed-file metadata read back after immutable scientific release')
        elif record:
            task.update(execution_status=record['status'],actual_output_paths=[record.get('directory',str(Path(output)/(str(stage)+'_latest_failure.json')))],block_reason=record.get('error'),stage=stage,fingerprint=record.get('fingerprint'))
            if record.get('directory'):
                completion=Path(record['directory'])/'completion.json';task['completion_manifest_path']=str(completion);task['completion_manifest_sha256']=digest(completion)
        elif tid=='T18' and index['stages'].get('figures',{}).get('directory'):
            task.update(execution_status='RUNNING',actual_output_paths=[index['stages']['figures']['directory']],scope='exports exist; actual visual inspection/verification pending')
    atomic_json(Path(output)/'TASK_STATUS.json',tasks)
    (Path(output)/'TASK_STATUS.tsv').write_text('id\ttitle\tevidence_status\texecution_status\n'+'\n'.join('\t'.join([t['id'],t['title'],t['evidence_status'],t['execution_status']]) for t in tasks['tasks'])+'\n')
