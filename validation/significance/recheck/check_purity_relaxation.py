"""Recheck four cached-kernel P1 profiles; no new quadrature or optimization.

By default use the bundled purity_inputs/ closure. An explicit --archive may
recover the same four members if that local closure is absent. Frozen results
are never overwritten: --output must name a new file.
"""
from pathlib import Path
import argparse
import hashlib
import io
import json
import sys
import tempfile
import zipfile

import numpy as np
from scipy.linalg import eigh

EXPECTED_ARCHIVE_SHA256='d311f1a1d84b15294379671d3b2a0346971d9aa2b4987f3dc6fc896b9d5b1222'
MEMBERS={'model_api_v1_1.py':'paper1_phase12_harmonic/model_api_v1_1.py',
         'reduced_core.py':'paper1_phase12_harmonic/reduced_core.py',
         'on_shell_N80.npz':'paper1_phase12_harmonic/raw/on_shell_N80.npz',
         'pv_mean_M192.npy':'paper1_phase12_harmonic/raw/pv_mean_M192.npy'}
MATRIX_TOL=1e-12
PSD_TOL=1e-9  # Numerical consistency only, as in original on-shell PSD test.
FLUX_TOL=1e-9
PROFILE_SUM_TOL=1e-14


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def inputs_at(path,archive):
    provenance=path/'INPUT_PROVENANCE.json'
    if provenance.is_file():
        metadata=json.loads(provenance.read_text())
        if metadata['source_archive_sha256']!=EXPECTED_ARCHIVE_SHA256:
            raise ValueError('Wrong source archive binding')
        if {x['path'] for x in metadata['files']}!=set(MEMBERS):
            raise ValueError('Wrong small input closure')
        payloads={}
        for item in metadata['files']:
            if item['archive_member']!=MEMBERS[item['path']]:
                raise ValueError('Wrong source archive member')
            p=path/item['path']
            if not p.is_file():raise ValueError('Missing input: '+item['path'])
            payload=p.read_bytes()
            if len(payload)!=item['bytes'] or digest(payload)!=item['sha256']:
                raise ValueError('Changed input: '+item['path'])
            payloads[item['path']]=payload
        return payloads,metadata,'bundled_local_inputs'
    if archive is None:
        raise ValueError('Local purity_inputs closure unavailable; provide explicit --archive only for exact original fallback')
    if digest(archive.read_bytes())!=EXPECTED_ARCHIVE_SHA256:
        raise ValueError('Fallback archive hash mismatch')
    with zipfile.ZipFile(archive) as z:
        payloads={name:z.read(member) for name,member in MEMBERS.items()}
    return payloads,{'source_archive_sha256':EXPECTED_ARCHIVE_SHA256,
       'files':[{'path':name,'archive_member':MEMBERS[name],'bytes':len(v),'sha256':digest(v)} for name,v in payloads.items()]},'explicit_archive_fallback'


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,default=here/'purity_inputs')
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Output exists; choose a fresh report filename')
    payloads,metadata,origin=inputs_at(args.inputs,args.archive)
    with np.load(io.BytesIO(payloads['on_shell_N80.npz']),allow_pickle=False) as data:
        d={k:data[k].copy() for k in data.files}
    pv=np.load(io.BytesIO(payloads['pv_mean_M192.npy']),allow_pickle=False)
    sl=np.s_[:12,:12]
    B0=d['Im_branch_0'][sl];B1=d['Im_branch_1'][sl]
    BL=d['Im_branch_2'][sl];BT=B0+B1;B=BT+BL
    total=d['Im_total'][sl];G=pv[sl]+1j*B
    matrices={'lower_pair_qT':BT,'upper_branch_qL':BL,'total':total,'PV':pv[sl]}
    matrix_checks={}
    for name,M in matrices.items():
        residual=float(np.linalg.norm(M-M.conj().T))
        assert residual<MATRIX_TOL,(name,residual)
        entry={'hermiticity_residual':residual}
        if name!='PV':
            minimum=float(np.linalg.eigvalsh((M+M.conj().T)/2).min())
            assert minimum>-PSD_TOL,(name,minimum)
            entry['minimum_eigenvalue']=minimum
        matrix_checks[name]=entry
    branch_sum_error=float(np.linalg.norm(B-total))
    assert branch_sum_error<MATRIX_TOL,branch_sum_error
    geometric_difference=float(np.linalg.norm(d['Im_T'][sl]-BT))
    assert geometric_difference>1e-7,'Naming adverse control not discriminating'

    def relaxation(W):
        T=W.conj().T@BT@W;L=W.conj().T@BL@W;S=T+L
        T=(T+T.conj().T)/2;L=(L+L.conj().T)/2;S=(S+S.conj().T)/2
        e,V=eigh(S);keep=e>1e-11*e.max();H=V[:,keep]/np.sqrt(e[keep])[None,:]
        rates=eigh(H.conj().T@T@H,eigvals_only=True)
        assert rates.min()>-PSD_TOL and rates.max()<1+PSD_TOL
        return {'radiating_rank':int(keep.sum()),'smallest_total_eigenvalue':float(e.min()),
                'qL_eigen_min':float(np.linalg.eigvalsh(L).min()),
                'maximum_relaxed_qT_share':float(rates.max()),'minimum_relaxed_qT_share':float(rates.min())}

    results={'status':'PASS','input_origin':origin,'input_provenance':metadata,'omega':1.2,'nseg':6,
      'cached_member_sha256':{k:digest(v) for k,v in payloads.items()},
      'numerical_controls':{'branch_sum_versus_Im_total_error':branch_sum_error,
         'geometric_Im_T_versus_lower_pair_norm':geometric_difference,'matrices':matrix_checks,
         'tolerances':{'matrix_identity_and_Hermiticity_absolute':MATRIX_TOL,'numerical_PSD_and_flux_absolute':PSD_TOL,'profile_sum_absolute':PROFILE_SUM_TOL}},
      'all_12_source_vectors':relaxation(np.eye(12)),'profiles':[]}
    profiles={'reference':np.array([.1,.2,.2,.2,.2,.1]),'uniform':np.ones(6)/6,
      'left_heavy':np.array([.4,.2,.15,.1,.1,.05]),'sparse_endpoints':np.array([.5,0,0,0,0,.5])}
    # Use exact original producer modules in a disposable location; no bytecode or
    # generated reports are written into bundled immutable input files.
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-purity-originals-') as td:
        stage=Path(td)
        for name in ['model_api_v1_1.py','reduced_core.py']:(stage/name).write_bytes(payloads[name])
        sys.path.insert(0,str(stage))
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices,scatter_from_reduced_green
        p=ModelParameters()
        for name,g in profiles.items():
            assert g.shape==(6,) and np.all(g>=0)
            total_strength=float(g.sum());assert abs(total_strength-1)<PROFILE_SUM_TOL
            U,V,_=reduced_static_matrices(g,p);W=np.linalg.solve(np.eye(12)+U@G,V)
            s,x=scatter_from_reduced_green(1.2,g,p,G,return_internal=True)
            f=x['t'];factor=2*x['gamma_in']
            qT=float(factor*np.vdot(f,BT@f).real);qL=float(factor*np.vdot(f,BL@f).real)
            bulk=float(factor*np.vdot(f,total@f).real)
            component_error=abs(qT+qL-bulk)
            flux_error=abs(s.Rq+s.Tq+qT+qL-1)
            assert component_error<MATRIX_TOL and abs(bulk-s.Tbulk)<MATRIX_TOL
            assert flux_error<FLUX_TOL and s.closure_error<FLUX_TOL
            assert s.min_gamma_bulk_eigenvalue>-PSD_TOL
            assert min(qT,qL,s.Rq,s.Tq)>-PSD_TOL
            results['profiles'].append({'name':name,'g':g.tolist(),'total_strength':total_strength,
              'arbitrary_line_coordinate_relaxation':relaxation(W),'actual_left_incident_qT':qT,
              'actual_left_incident_qL':qL,'actual_qT_share':qT/(qT+qL),'line_reflection':s.Rq,
              'line_transmission':s.Tq,'closure':s.closure_error,'independent_four_channel_closure':flux_error,
              'qT_plus_qL_versus_total_error':component_error,'minimum_loaded_escape_eigenvalue':s.min_gamma_bulk_eigenvalue})
    results['limits']=['Recheck of the same four exploratory historical N80 profiles, not new parameter points or current P1 coefficient updates.',
      'Numerical Hermiticity/PSD/flux checks are internal consistency, not exact infinite-lattice enclosures.',
      'qT is lower two BRANCH sum; geometric Im_T is deliberately not used as qT.',
      'Relaxed maxima are not proven physically attainable by a positive passive profile or scalar incident state.',
      'No fresh infinite-lattice quadrature, production trajectory, optimization or field-wide novelty claim.']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(results,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'status':'PASS','report':str(args.output),'profiles_checked':len(results['profiles']),
      'branch_sum_error':branch_sum_error,'maximum_four_channel_closure':max(x['independent_four_channel_closure'] for x in results['profiles'])},indent=2))


if __name__=='__main__':
    main()
