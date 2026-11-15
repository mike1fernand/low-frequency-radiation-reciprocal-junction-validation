"""Independently bind the point center to replayed R4 exact targets.

This verifies target identity and safe rational caps. It consumes the separately
executed full R4 saved-evidence replay; it does not rerun its global quadrature.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json
import numpy as np

BASE=Path('/workspace/scratch/044aa7b17d80/p1_736_validation_2026_11_08')
OUT=Path(__file__).resolve().parent
KERNEL=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/REFERENCE_KERNEL.npz')
BUDGET=BASE/'r4_certificate/analysis/full_domain_r4/budget'
targets_path=BUDGET/'TARGETS.json'
ledger_path=BUDGET/'FINAL_CAMPAIGN_LEDGER.json'
replay_path=BASE/'baseline_audit/r4_replay_out/VERIFICATION.json'
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
targets=json.loads(targets_path.read_text())
ledger=json.loads(ledger_path.read_text())
replay=json.loads(replay_path.read_text())
assert replay['status']=='passed'
assert replay['exact_lattice_floor_certified'] is True
assert replay['physical_replay_performed'] is True
assert replay['saved_ledger_reproduced_exactly'] is True
assert all(x['status']=='passed' and x['returncode']==0 for x in replay['stages'])
assert replay['ledger_sha256']==H(ledger_path)
assert replay['bridge']['bindings']['targets_sha256']==H(targets_path)
assert H(KERNEL)=='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
with np.load(KERNEL,allow_pickle=False) as z:
    G=z['G0'][:12,:12]
    Q=z['Im_qT'][:12,:12]
    B0=z['Im_branch_0'][:12,:12]
    B1=z['Im_branch_1'][:12,:12]
for i in range(12):
    for j in range(12):
        assert F(float(G[i,j].real))==F(targets['targets']['g'][i][j][0])
        assert F(float(G[i,j].imag))==F(targets['targets']['g'][i][j][1])
diff=F(0);target_is_exact_branchsum=True
for i in range(12):
    for j in range(12):
        qt=F(targets['targets']['Q'][i][j][0])
        assert F(targets['targets']['Q'][i][j][1])==0
        target_is_exact_branchsum &= qt==F(float(B0[i,j]))+F(float(B1[i,j]))
        diff=max(diff,abs(F(float(Q[i,j]))-qt))
eG=F(ledger['results']['g']['operator_budget']['scaled_operator_upper'])
eQ=F(ledger['results']['Q']['operator_budget']['scaled_operator_upper'])
centering_allowance=F(169,100)*12*diff
assert eG<F(1,250)
assert eQ+centering_allowance<F(1,500)
result={
    'status':'PASS','date':'2026-11-08',
    'scope':'Exact target/center binding and upward norm-cap arithmetic, consuming the separately passed saved-evidence R4 replay. No fresh global contour integration.',
    'hashes':{'kernel':H(KERNEL),'targets':H(targets_path),
              'ledger':H(ledger_path),'saved_evidence_replay':H(replay_path)},
    'G_exact_center_matches_targets':True,
    'Q_target_matches_exact_sum_of_branch_dyads':target_is_exact_branchsum,
    'Q_float_sum_to_certified_target_max_entry':str(diff),
    'scaled_Q_centering_operator_allowance':str(centering_allowance),
    'eG_replayed_exact_rational':str(eG),
    'eQ_replayed_exact_rational':str(eQ),
    'safe_caps':{'eG':'1/250','eQ':'1/500'},
    'safe_caps_cover_center_shift':True,
    'saved_evidence_status':{'all_stages_passed':True,'physical_controls_and_bridge_passed':True,
                             'saved_ledger_reproduced_exactly':True},
}
(OUT/'KERNEL_CENTER_BINDING_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['eG_replayed_exact_rational','eQ_replayed_exact_rational']},indent=2))
