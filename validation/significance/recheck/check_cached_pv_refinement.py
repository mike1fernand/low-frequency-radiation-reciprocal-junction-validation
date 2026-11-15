"""Two previously seen profiles against a same-archive alternate PV cache only."""
from pathlib import Path
import argparse,hashlib,io,json,sys,tempfile
import numpy as np
import check_purity_relaxation as cached


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',type=Path,default=here/'purity_inputs')
    parser.add_argument('--refinement-inputs',type=Path,default=here/'purity_refinement_inputs');parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Choose a fresh output filename')
    payloads,provenance,origin=cached.inputs_at(args.inputs,None)
    context=json.loads((args.refinement_inputs/'INPUT_PROVENANCE.json').read_text())
    assert context['source_archive_sha256']==cached.EXPECTED_ARCHIVE_SHA256
    for item in context['files']:
        data=(args.refinement_inputs/item['path']).read_bytes()
        assert len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256']
    pv=np.load(args.refinement_inputs/'pv_direct_eta_extrapolated.npy',allow_pickle=False)[:12,:12]
    with np.load(io.BytesIO(payloads['on_shell_N80.npz']),allow_pickle=False) as d:
        BT=d['Im_branch_0'][:12,:12]+d['Im_branch_1'][:12,:12];BL=d['Im_branch_2'][:12,:12];B=d['Im_total'][:12,:12]
    assert np.linalg.norm(BT+BL-B)<1e-12
    old=json.loads((here/'PURITY_RELAXATION_RESULT.json').read_text())
    old_rows={x['name']:x for x in old['profiles']}
    out={'status':'PASS','input_provenance':provenance,'alternate_PV_provenance':context,'profiles':[]}
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-pv-control-originals-') as td:
        p=Path(td)
        for name in ['model_api_v1_1.py','reduced_core.py']:(p/name).write_bytes(payloads[name])
        sys.path.insert(0,str(p))
        from model_api_v1_1 import ModelParameters
        from reduced_core import scatter_from_reduced_green
        for name in ['reference','sparse_endpoints']:
            g=np.array(old_rows[name]['g']);assert abs(g.sum()-1)<1e-14
            r,x=scatter_from_reduced_green(1.2,g,ModelParameters(),pv+1j*B,return_internal=True)
            f=x['t'];qT=float(2*x['gamma_in']*np.vdot(f,BT@f).real);qL=float(2*x['gamma_in']*np.vdot(f,BL@f).real)
            flux=abs(r.Rq+r.Tq+qT+qL-1);assert flux<1e-9
            assert abs(qT+qL-r.Tbulk)<1e-12
            out['profiles'].append({'name':name,'g':g.tolist(),'incident_qT_alternate_PV':qT,'incident_qL_alternate_PV':qL,
               'qT_share_alternate_PV':qT/(qT+qL),'original_PV_incident_qT':old_rows[name]['actual_left_incident_qT'],
               'qT_absolute_PV_difference':abs(qT-old_rows[name]['actual_left_incident_qT']),'four_channel_closure':flux})
    out['endpoint_minus_reference_qT_alternate_PV']=out['profiles'][1]['incident_qT_alternate_PV']-out['profiles'][0]['incident_qT_alternate_PV']
    out['limits']=['Exactly two rechecks of previously seen profiles using an existing same-archive alternate PV estimate; no new quadrature.',
      'The on-shell N80 branch matrices are SHARED, so this is not an independent complete method or enclosure.',
      'No beta endpoint-imbalance state was refined, optimized or added.',
      'The existing historical direct-eta extrapolation metadata is retained; its difference is sensitivity, not a statistical confidence interval.']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
