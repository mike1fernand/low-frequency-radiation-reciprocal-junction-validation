"""Review bindings and strict inequalities in the fresh original bridge output."""
import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ctx.prec=160
    directory=ROOT/'bridge_replay'
    summary=json.loads((directory/'EXACT_LATTICE_BRIDGE.json').read_text())
    proof=json.loads((directory/'finite_proof/finite_kernel_certificate.json').read_text())
    gate=json.loads((directory/'GLOBAL_GATE.json').read_text())
    readiness=json.loads((ROOT/'budget/BRIDGE_READINESS.json').read_text())
    assert sha(directory/'GLOBAL_GATE.json')==summary['bindings']['global_gate_sha256']==sha(ROOT/'budget/FINAL_CAMPAIGN_LEDGER.json')
    assert sha(directory/'finite_proof/finite_kernel_certificate.json')==summary['bindings']['fresh_original_proof_sha256']
    assert sha(ROOT/'budget/bridge_inputs/finite_kernel_certificate.py')==summary['bindings']['original_producer_sha256']==readiness['sha256']['finite_kernel_certificate.py']
    assert sha(ROOT/'budget/TARGETS.json')==summary['bindings']['targets_sha256']
    baseline=(ROOT/json.loads((ROOT/'CAMPAIGN.json').read_text())['baseline']).resolve()
    assert sha(baseline/'math/COMPLETE_FALLBACK_PROOF.md')==summary['bindings']['exact_Q_psd_proof_sha256']
    assert {tuple(v['signs']) for v in proof['vertices']}==set(product((-1,1),repeat=6)) and len(proof['vertices'])==64
    assert proof['status']=='finite_cached_kernel_certificate_passed'
    bridge=proof['hypothetical_kernel_error_conclusion']
    assert bridge['eG']=='1/50' and bridge['eQ']=='1/200' and bridge['rational_floor']=='1037/10000'
    assert arb(bridge['floor_enclosure'])>arb(1037)/10000
    assert arb(bridge['full_inverse_upper'])<1 and arb(bridge['core_inverse_upper'])<1
    assert all(F(x['fallback_parameter_area'])==0 for x in gate['results'].values())
    for kind,threshold in [('g',F(1,50)),('Q',F(1,200))]:
        assert F(gate['results'][kind]['operator_budget']['direct_envelope']['scaled_operator_upper'])<threshold
    result={'status':'passed','checked_vertices':64,
            'checked_positive_LDL_pivots':proof['checked_positive_LDL_pivots'],
            'bridge_floor':bridge['floor_enclosure'],
            'full_inverse_upper':bridge['full_inverse_upper'],
            'core_inverse_upper':bridge['core_inverse_upper'],
            'summary_sha256':sha(directory/'EXACT_LATTICE_BRIDGE.json'),
            'proof_sha256':sha(directory/'finite_proof/finite_kernel_certificate.json'),
            'ledger_sha256':sha(directory/'GLOBAL_GATE.json'),
            'scope':'Reference energy 36/25, specified source, all six complete-motif coupling factors in [19/20,21/20].',
            'method':'Independently checked all bindings, exact set of 64 vertex signs, positive floor margin and both inverse conditions in the fresh original proof. Original proof was not rerun a second time by this reviewer.'}
    Path(__file__).with_name('BRIDGE_REVIEW_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
