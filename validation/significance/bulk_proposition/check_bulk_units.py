"""Check the newly restored physical-unit/domain statement against the engine.

Nine small static cases, not a repeat of the old 35+320 checks, a trajectory,
quadrature, global proof, or validation of generated power.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json
import numpy as np


def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,default=here.parent/'recheck'/'purity_inputs')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():ap.error('Choose a new output filename')
    source=args.inputs/'model_api_v1_1.py'
    expected='aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980'
    assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
    # No bytecode/files are written into the preserved engine directory.
    namespace={'__name__':'bulk_units_archived_model'}
    import types,sys
    module=types.ModuleType(namespace['__name__']);sys.modules[module.__name__]=module
    exec(compile(source.read_text(),str(source),'exec'),module.__dict__)
    bonds=np.asarray([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),
                      (1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
    cases=[];wrong_unit_difference=0.
    for a,m,K in [(1.,1.,1.),(.25,2.,3.),(2.,.5,.7)]:
        p=module.ModelParameters(a=a,m_bulk=m,k_axis=K,k_diag=K)
        for rho,direction in [(.2,(1.,0.,0.)),(1.1,(1.,2.,3.)),(1.99,(1.,1.,0.))]:
            n=np.asarray(direction);n/=np.linalg.norm(n);xi=rho*n;kphysical=xi/a
            physical_bonds=a*bonds;unit=bonds/np.linalg.norm(bonds,axis=1)[:,None]
            weights=(2*K/m)*(1-np.cos(physical_bonds@kphysical))
            direct=np.einsum('b,bi,bj->ij',weights,unit,unit)
            engine=module.bulk_dynamical_matrix(xi,p)
            scale=max(1.,np.linalg.norm(direct));err=np.linalg.norm(engine-direct)/scale
            assert err<2e-14,err
            Dhat=(m/K)*engine;D2=rho*rho*np.eye(3)+2*np.outer(xi,xi)
            E=D2-Dhat;ev=np.linalg.eigvalsh(E)
            assert ev.min()>-1e-12 and ev.max()<=rho**4/2+1e-12
            values,vectors=np.linalg.eigh(engine);P=vectors[:,:2]@vectors[:,:2].T
            lower_gap=(K/m)*(2*rho*rho-rho**4/2)
            gap=values[2]-values[1]
            assert lower_gap>0 and gap>=lower_gap-1e-12
            amplitude=float(np.linalg.norm(P@n));leakage=float(n@P@n)
            assert amplitude<=rho*rho/4+1e-12
            assert abs(amplitude*amplitude-leakage)<1e-13
            assert leakage<=rho**4/16+1e-12
            if a!=1.:
                wrong=module.bulk_dynamical_matrix(kphysical,p)
                wrong_unit_difference=max(wrong_unit_difference,float(np.linalg.norm(wrong-direct)))
            cases.append({'a':a,'mass':m,'stiffness':K,'rho':rho,
                          'relative_engine_symbol_error':float(err),'physical_gap':float(gap),
                          'sufficient_gap_lower_bound':lower_gap,'amplitude_leakage':amplitude,
                          'squared_projector_leakage':leakage})
    assert wrong_unit_difference>1e-3,'Adverse unit-convention control must discriminate'
    out={'status':'PASS','scope':'Nine targeted static physical-unit calibrations. Not a global proof, reference-shell enclosure, generated-power result or production simulation.',
         'cases':cases,'wrong_units_adverse_max_difference':wrong_unit_difference,
         'engine_sha256':expected,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'existing_35_plus_320_cases_rerun':False,'new_time_domain_simulations':0}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':'PASS','cases':len(cases),
                      'maximum_symbol_error':max(x['relative_engine_symbol_error'] for x in cases),
                      'wrong_units_adverse_difference':wrong_unit_difference}))


if __name__=='__main__':main()
