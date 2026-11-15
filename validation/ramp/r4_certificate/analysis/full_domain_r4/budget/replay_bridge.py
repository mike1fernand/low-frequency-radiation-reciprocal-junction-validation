"""Run the unchanged selected-matrix bridge after actual global budgets pass.

The full original proof is cheap and not modular; no replacement derivation
or rounded scalar substitution is introduced here. Failure of any prerequisite
leaves the existing finite-matrix claim unchanged.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import subprocess
import sys

import campaign_ledger
import matrix_budget as mb


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();output=a.output_dir.resolve()
    if output.exists():raise ValueError('Fresh output directory required')
    package=Path(__file__).resolve().parent
    inputs=package/'bridge_inputs'
    readiness=json.loads((package/'BRIDGE_READINESS.json').read_text())
    for relative,expected in readiness['sha256'].items():
        if sha(inputs/relative)!=expected:raise ValueError('Original bridge input/source binding mismatch: '+relative)
    # Recompute the gate from audited checkpoints. Do not accept cached booleans.
    ledger=campaign_ledger.assemble(root,root/'CAMPAIGN.json')
    if not ledger['both_arithmetic_budgets_pass']:
        raise ValueError('Both unchanged actual global operator budgets must pass before bridge replay')
    if any(F(x['fallback_parameter_area']) for x in ledger['results'].values()):
        raise ValueError('Complete refined campaign required for this bridge promotion')
    output.mkdir(parents=True)
    gate_path=output/'GLOBAL_GATE.json'
    gate_path.write_text(json.dumps(mb.json_exact(ledger),indent=2)+'\n')
    cmd=[sys.executable,'-B',str(inputs/'finite_kernel_certificate.py'),
         '--repo-root',str(inputs/'staged_repo'),'--output-dir',str(output/'finite_proof'),
         '--precision','160']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=45)
    (output/'finite_proof.stdout.txt').write_text(result.stdout)
    (output/'finite_proof.stderr.txt').write_text(result.stderr)
    if result.returncode!=0:raise RuntimeError('Original bridge replay did not pass; see saved process outputs')
    proof_path=output/'finite_proof/finite_kernel_certificate.json'
    proof=json.loads(proof_path.read_text())
    bridge=proof['hypothetical_kernel_error_conclusion']
    if (proof['status']!='finite_cached_kernel_certificate_passed' or
            proof['certificate_revision']!='basis_invariant_n192_sobol23_m256_v1' or
            proof['finite_rational_floor']!='27/250' or
            bridge['eG']!='1/50' or bridge['eQ']!='1/200' or
            bridge['rational_floor']!='1037/10000' or len(proof['vertices'])!=64):
        raise ValueError('Fresh original proof did not match the stated bridge contract')
    targets=json.loads((package/'TARGETS.json').read_text())
    needed={key:readiness['sha256']['staged_repo/'+key] for key in proof['source_sha256']}
    if proof['source_sha256']!=needed:raise ValueError('Fresh bridge targets changed')
    baseline=(root/json.loads((root/'CAMPAIGN.json').read_text())['baseline']).resolve()
    psd=baseline/'math/COMPLETE_FALLBACK_PROOF.md'
    summary={'status':'exact_reference_lattice_floor_certified_from_complete_enclosures_and_original_bridge',
             'exact_lattice_floor_certified':True,'rational_floor':'1037/10000',
             'claim':'At the specified reference energy E=36/25 and the specified source, T_qT>1037/10000 throughout the same six independent complete-motif coupling-factor box [19/20,21/20].',
             'scope_limits':'No packet-averaged guarantee, no frequency-band theorem, and no experimental/model-validity certification.',
             'psd_premise':'Exact-Q PSD follows from the complete lower-pair spectral-residue quadratic form in the bound inherited proof; no midpoint PSD inference is used.',
             'bindings':{'global_gate_sha256':sha(gate_path),'fresh_original_proof_sha256':sha(proof_path),
                         'original_producer_sha256':sha(inputs/'finite_kernel_certificate.py'),
                         'exact_Q_psd_proof_sha256':sha(psd),'targets_sha256':sha(package/'TARGETS.json')},
             'fresh_original_process_exit':result.returncode,'original_conditional_floor_enclosure':bridge['floor_enclosure']}
    (output/'EXACT_LATTICE_BRIDGE.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
