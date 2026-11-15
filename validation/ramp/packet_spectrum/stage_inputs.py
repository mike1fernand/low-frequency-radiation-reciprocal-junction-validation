"""Stage only the exact fourteen required packet nodes from the recovered R8 archive."""
from pathlib import Path
import hashlib,json,shutil

HERE=Path(__file__).resolve().parent
TOP=HERE.parent
BASE=TOP/'r8_selected/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics/analysis/resolution_execution'
WORK=TOP/'r8_selected/paper1-validation-integrated-r8/legacy_rc4/validation/remaining_resolution/packet/worktree/analysis/resolution_execution'
INPUT=HERE/'inputs';INPUT.mkdir(parents=True,exist_ok=True)
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,x:p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
models=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs')
copied={}
for n in ['model_api_v1_1.py','reduced_core.py']:
    shutil.copyfile(models/n,INPUT/n);copied[n]=H(INPUT/n)
for n in ['packet_definition.py','FREQUENCY_CONTRACT.json']:
    shutil.copyfile(BASE/'packet'/n,INPUT/n);copied[n]=H(INPUT/n)
shutil.copyfile(BASE/'packet/inputs/preparation.json',INPUT/'preparation.json');copied['preparation.json']=H(INPUT/'preparation.json')
contract=json.loads((INPUT/'FREQUENCY_CONTRACT.json').read_text())
assert copied['preparation.json']==contract['input_sha256']
assert contract['central_node_count']==14
source_by_hash={}
for root in [BASE,WORK]:
    for p in root.rglob('*.py'):source_by_hash.setdefault(H(p),p)
families=[]
for label,method,n,root in [
 ('surface_N128','surface',128,BASE/'surface/selected_nodes/N128'),
 ('surface_N160','surface',160,BASE/'surface/selected_nodes/N160'),
 ('tetra_N96','tetrahedron',96,BASE/'tetra/packet_N96'),
 ('tetra_N128','tetrahedron',128,BASE/'tetra/packet_N128'),
 ('tetra_N160','tetrahedron',160,WORK/'tetra/packet_N160_complete')]:
    dest=INPUT/label;dest.mkdir(exist_ok=True)
    records=[]
    for p in root.glob('omega*.json'):
        d=json.loads(p.read_text())
        if 'omega' in d and p.with_suffix('.npz').is_file():records.append((p,d))
    family={'label':label,'method':method,'grid':n,'source_directory':str(root),'nodes':[]}
    if method=='tetrahedron':
        cfg=json.loads((root/'CONFIG.json').read_text());binding=json.loads((root/'PACKET_FREQUENCY_INPUT.json').read_text())
        assert cfg['config']['nseg']==6 and cfg['config']['profile']==[.1,.2,.2,.2,.2,.1]
        assert cfg['config']['grid']==n and cfg['config']['shift']==[0,0,0]
        assert binding['packet_contract_sha256']==copied['FREQUENCY_CONTRACT.json']
        assert binding['central_node_ids']==[r['node_id'] for r in contract['central_nodes']]
        source_checks=[]
        for rel,h in cfg['config']['source_sha256'].items():
            actual=root.parents[3]/rel
            assert actual.is_file() and H(actual)==h,(label,rel)
            source_checks.append({'path':rel,'sha256':h,'matches':True})
        family['bulk_regulator']=cfg['config']['regulator'];family['source_checks']=source_checks
        for nmeta in ['CONFIG.json','PACKET_FREQUENCY_INPUT.json','RESULT.json','LAST_EXECUTION.json']:
            shutil.copyfile(root/nmeta,dest/nmeta)
    for node in contract['central_nodes']:
        matches=[(p,d) for p,d in records if abs(d['omega']-node['omega'])<2e-13]
        assert len(matches)==1,(label,node,matches)
        p,meta=matches[0];kp=p.with_suffix('.npz')
        assert H(kp)==meta['kernel_sha256']
        assert meta['grid']==n
        if method=='surface':
            assert meta['support_segments']==6 and meta['profile']==[.1,.2,.2,.2,.2,.1]
            assert meta['frequency_binding']['node_id']==node['node_id']
            assert meta['frequency_binding']['packet_contract_sha256']==copied['FREQUENCY_CONTRACT.json']
            assert meta['producer_sha256'] in source_by_hash
            family.setdefault('selected_producer_checks',{})[meta['producer_sha256']]=str(source_by_hash[meta['producer_sha256']])
        local_npz=dest/(node['node_id']+'.npz');local_meta=dest/(node['node_id']+'.json')
        shutil.copyfile(kp,local_npz);shutil.copyfile(p,local_meta)
        family['nodes'].append({'node_id':node['node_id'],'omega':meta['omega'],
             'kernel':str(local_npz.relative_to(HERE)),'metadata':str(local_meta.relative_to(HERE)),
             'kernel_sha256':H(kp),'metadata_sha256':H(p),
             'archive_relative_kernel':str(kp.relative_to(TOP/'r8_selected')),
             'archive_relative_metadata':str(p.relative_to(TOP/'r8_selected'))})
    families.append(family)
manifest={'status':'PASS','scope':'Exact fourteen in-band positive-index DFT nodes for five already-computed bulk kernel families. Original reference and new ramp will be solved separately; no new Green quadrature.',
 'source_archive_sha256':'71259248c3e0c5faeb4e508962ae6d25693df55b1ebf642fc3b88ddd26e80b9e',
 'source_archive_library_id':'libfile_a5d9914430088191891ba57eddc56542','copied_sources':copied,'families':families}
J(HERE/'INPUT_BINDING.json',manifest)
print(json.dumps({'status':'PASS','families':len(families),'kernels':sum(len(x['nodes']) for x in families)},indent=2))
