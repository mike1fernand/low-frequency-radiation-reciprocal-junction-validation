"""Small exact normalization check; no wave simulation or new source quadrature."""
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,math
# Pin=(Z/2)omega^2 |U|^2, |F|^2=4 Z^2 omega^2 |U|^2,
# P=(omega/2)|F|^2 Im(chi), Im(chi)=omega*b.
coefficient=Q(1,2)*Q(4)/Q(1,2)
assert coefficient==4
assert coefficient*Q(1,6)==Q(2,3)
assert coefficient*Q(1,12)==Q(1,3)
# Published numbers are rounded; this is a conversion consistency check only,
# not a precision recalculation of their half-space integral.
p=.694
r2=(3-4*p*p)/(1-p*p)
converted=1.91/r2
assert abs(converted-.923)<.001
out={'status':'PASS','scope':'Exact time-average/force normalization; rounded published preprint coefficient conversion only. No simulation, full-space/half-space equivalence, independent half-space integral or novelty proof.',
'exact_T_over_Z_b_omega_squared':str(coefficient),
'exact_qT_multiplier_without_pi_rho_cubedspeed':'2/3',
'exact_qL_multiplier_without_pi_rho_cubedspeed':'1/3',
'CG_preprint_p':p,'CG_c0_squared_over_ct_squared':r2,
'CG_rounded_1_91_bk_squared_converted_to_bomega_over_ct_squared':converted,
'CG_printed_rounded_value':.923,
'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
