"""Verify compact payloads and authenticated proof-chain input bindings.

This checks saved bytes and statuses; it does not replace re-executing the
rational response verifier or the complete baseline certificate verifier.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json

ROOT=Path(__file__).resolve().parent
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((ROOT/'PACKAGE_SHA256SUMS.json').read_text())
for name,expected in manifest.items():
    p=ROOT/name
    assert p.is_file() and H(p)==expected,'Payload mismatch: '+name
chain=json.loads((ROOT/'tolerance/VALIDATED_PROOF_CHAIN.json').read_text())
for name,binding in chain['files'].items():
    p=ROOT/binding['path']
    assert p.is_file() and H(p)==binding['sha256'],'Proof binding mismatch: '+name
point=chain['nominal_infinite_lattice_interval']
box=chain['whole_six_motif_plus_minus5_percent_infinite_lattice_floor']
assert Fraction(point['lower_rational'])>Fraction(732,1000)
assert Fraction(point['upper_rational'])<Fraction(740,1000)
assert Fraction(box['lower_rational'])>Fraction(69,100)
assert chain['upstream_status']['saved_ledger_reproduced_exactly']
assert chain['upstream_status']['combined_replay_with_physical_and_bridge_status_certified']
print(json.dumps({'status':'PASS_SAVED_BYTES_AND_PROOF_BINDINGS',
                  'payloads':len(manifest),'proof_chain_bindings':len(chain['files']),
                  'full_baseline_source_present':(ROOT/'r4_certificate/MANIFEST.json').exists(),
                  'reexecution_note':'Run the separately documented response and full-baseline verifiers for independent proof re-execution.'},indent=2))
