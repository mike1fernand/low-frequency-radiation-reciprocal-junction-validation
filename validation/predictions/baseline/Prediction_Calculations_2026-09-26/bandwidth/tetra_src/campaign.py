#!/usr/bin/env python3
"""Resume checked cube-range checkpoints; frequency rows share frozen geometry."""
from pathlib import Path
from dataclasses import asdict
import argparse,json,hashlib,os,time
import numpy as np
import multifrequency_tetra as m

def atomic_json(path,data):
    tmp=Path(str(path)+'.pending');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)

def load_part(directory,start,stop,config_hash,nf,nseg=6):
    base=directory/f'part_{start:09d}_{stop:09d}';meta=json.loads(base.with_suffix('.json').read_text());fp=base.with_suffix('.npz')
    if meta['config_hash']!=config_hash or meta['start']!=start or meta['stop']!=stop or m.sha(fp)!=meta['npz_sha256']:raise ValueError('checkpoint identity/range/hash mismatch')
    with np.load(fp,allow_pickle=False) as x:
        sums=x['sums'];mins=x['minimum_scalar_imaginary'];neg=x['negative_scalar_count'];maxs=x['maximum_abs_scalar']
    if sums.shape!=(nf,3,2,nseg,2,2) or mins.shape!=(nf,) or neg.shape!=(nf,) or maxs.shape!=(nf,):raise ValueError('checkpoint shape mismatch')
    if not all(np.isfinite(x).all() for x in [sums,mins,neg,maxs]) or np.any(mins<0) or np.any(neg!=0):raise ValueError('checkpoint nonfinite/negative spectral factors')
    return dict(meta,sums=sums,minimum_scalar_imaginary=mins,negative_scalar_count=neg,maximum_abs_scalar=maxs)

def _campaign(directory,grid,omegas,shift=(0,0,0),regulator=1e-10,range_cubes=65536,cube_chunk=1024,max_new_parts=None,merge_only=False):
    if int(range_cubes)!=range_cubes or range_cubes<1 or int(cube_chunk)!=cube_chunk or cube_chunk<1:raise ValueError('positive integer cube counts required')
    if max_new_parts is not None and (int(max_new_parts)!=max_new_parts or max_new_parts<0):raise ValueError('nonnegative integer checkpoint limit required')
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    ev,geometry=m.prepare_geometry(directory/'geometry',grid,shift)
    source=m.source_hashes()|{str(Path(__file__).relative_to(m.ROOT)):m.sha(__file__)}
    config={'grid':grid,'omega':[float(x) for x in omegas],'shift':[float(x) for x in shift],'regulator':regulator,'range_cubes':range_cubes,'cube_chunk':cube_chunk,'nseg':6,'profile':[.1,.2,.2,.2,.2,.1],'parameters':asdict(m.accepted.Parameters()),'geometry':geometry['settings'],'geometry_sha256':geometry['eigenvalues_sha256'],'source_sha256':source}
    digest=hashlib.sha256(m.canonical(config).encode()).hexdigest();cp=directory/'CONFIG.json'
    if cp.exists() and json.loads(cp.read_text())!={'config_hash':digest,'config':config}:raise ValueError('campaign configuration mismatch')
    if not cp.exists():atomic_json(cp,{'config_hash':digest,'config':config})
    ranges=[(s,min(s+range_cubes,grid**3)) for s in range(0,grid**3,range_cubes)];loaded=0;created=0;parts=[]
    for start,stop in ranges:
        base=directory/f'part_{start:09d}_{stop:09d}';mf=base.with_suffix('.json')
        if mf.exists():part=load_part(directory,start,stop,digest,len(omegas));loaded+=1
        else:
            if merge_only:raise FileNotFoundError('required checkpoint missing: '+str(mf))
            if max_new_parts is not None and created>=max_new_parts:
                progress={'status':'partial','config_hash':digest,'loaded':loaded,'created':created,'required':len(ranges),'next_start':start};atomic_json(directory/'LAST_EXECUTION.json',progress);return progress
            part=m.integrate_range(ev,omegas,grid,start,stop,shift=shift,regulator=regulator,cube_chunk=cube_chunk)
            fp=base.with_suffix('.npz');tmp=Path(str(base)+'.pending.npz');np.savez_compressed(tmp,**{k:part[k] for k in ['sums','minimum_scalar_imaginary','negative_scalar_count','maximum_abs_scalar']});os.replace(tmp,fp)
            meta={k:part[k] for k in ['centroid_clusters','seconds','start','stop']}|{'config_hash':digest,'npz_sha256':m.sha(fp)};atomic_json(mf,meta)
            part=load_part(directory,start,stop,digest,len(omegas));created+=1
            print(json.dumps({'checkpoint':[start,stop],'seconds':part['seconds'],'config_hash':digest}),flush=True)
        parts.append(part)
    # Exactly-once, ordered range coverage was enforced by the fixed range list;
    # every row is independently checked against the same frozen configuration.
    summed=np.zeros_like(parts[0]['sums']);mins=np.full(len(omegas),np.inf);clusters={}
    for part in parts:
        summed+=part['sums'];mins=np.minimum(mins,part['minimum_scalar_imaginary'])
        for k,v in part['centroid_clusters'].items():clusters[k]=clusters.get(k,0)+v
    rows=[]
    for i,omega in enumerate(omegas):
        diagnostics={'seconds':sum(x['seconds'] for x in parts),'minimum_scalar_imaginary':float(mins[i]),'centroid_clusters':clusters,'range_count':len(parts)}
        kernel=m.result_for_frequency(omega,grid,regulator,summed[i],diagnostics);scatter=m.accepted.independent_scatter(omega,config['profile'],m.accepted.Parameters(),kernel)
        if scatter.closure_amplitude>1e-8 or scatter.component_closure_modal>1e-10:raise ArithmeticError('scattering closure failed')
        base=directory/f'omega_{i:03d}';np.savez_compressed(base.with_suffix('.npz'),G0=kernel.G0,Im0=kernel.Im_branches[0],Im1=kernel.Im_branches[1],Im2=kernel.Im_branches[2],ImL=kernel.Im_longitudinal_projector,ImBL0=kernel.Im_branch_longitudinal_projectors[0],ImBL1=kernel.Im_branch_longitudinal_projectors[1],ImBL2=kernel.Im_branch_longitudinal_projectors[2])
        row={'omega':omega,'grid':grid,'shift':list(shift),'regulator':regulator,'diagnostics':kernel.diagnostics,'scattering':asdict(scatter),'kernel_sha256':m.sha(base.with_suffix('.npz'))};atomic_json(base.with_suffix('.json'),row);rows.append(row)
    result={'status':'complete_numerical_sensitivity_run_not_error_certificate','config_hash':digest,'geometry_seconds':geometry['seconds'],'shared_integration_seconds':sum(x['seconds'] for x in parts),'rows':rows}
    atomic_json(directory/'RESULT.json',result);atomic_json(directory/'LAST_EXECUTION.json',{'status':'complete','config_hash':digest,'loaded':loaded,'created':created,'required':len(ranges)})
    print(json.dumps({'status':result['status'],'grid':grid,'shift':list(shift),'seconds':result['shared_integration_seconds'],'qT':[r['scattering']['TT_modal'] for r in rows]}),flush=True);return result

def campaign(directory,grid,omegas,shift=(0,0,0),regulator=1e-10,range_cubes=65536,cube_chunk=1024,max_new_parts=None,merge_only=False):
    # Execution lifecycle is separate from preserved successful scientific data.
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    result_path=directory/'RESULT.json';last_successful=None
    if result_path.exists():
        try:last_successful=json.loads(result_path.read_text()).get('config_hash')
        except (ValueError,OSError):pass
    request={'grid':repr(grid),'omega':repr(omegas),'shift':repr(shift),'regulator':repr(regulator),'range_cubes':repr(range_cubes),'cube_chunk':repr(cube_chunk),'max_new_parts':repr(max_new_parts),'merge_only':bool(merge_only)}
    status_path=directory/'LAST_EXECUTION.json'
    atomic_json(status_path,{'status':'running','request':request,'previous_successful_result_config_hash':last_successful,'previous_result_is_from_earlier_execution':last_successful is not None})
    try:
        result=_campaign(directory,grid,omegas,shift,regulator,range_cubes,cube_chunk,max_new_parts,merge_only)
    except BaseException as error:
        atomic_json(status_path,{'status':'failed','request':request,'error_type':type(error).__name__,'error':str(error),'previous_successful_result_config_hash':last_successful,'previous_result_is_from_earlier_execution':last_successful is not None})
        raise
    state=json.loads(status_path.read_text());state['request']=request
    state['previous_result_is_from_earlier_execution']=state['status']!='complete' and last_successful is not None
    state['previous_successful_result_config_hash']=last_successful
    atomic_json(status_path,state)
    return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--directory',type=Path,required=True);a.add_argument('--grid',type=int,required=True);a.add_argument('--omega',nargs='+',type=float,required=True);a.add_argument('--shift',nargs=3,type=float,default=[0,0,0]);a.add_argument('--regulator',type=float,default=1e-10);a.add_argument('--range-cubes',type=int,default=65536);a.add_argument('--chunk',type=int,default=1024);a.add_argument('--max-new-parts',type=int);a.add_argument('--merge-only',action='store_true');q=a.parse_args()
    if q.range_cubes<1 or q.chunk<1:a.error('positive cube counts required')
    campaign(q.directory,q.grid,q.omega,q.shift,q.regulator,q.range_cubes,q.chunk,q.max_new_parts,q.merge_only)
