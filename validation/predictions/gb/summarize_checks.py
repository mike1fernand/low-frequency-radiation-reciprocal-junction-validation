import json,hashlib,platform
from pathlib import Path
import numpy as np
import p1solver as p
from independent_core_scattering import evaluate
R=Path(__file__).resolve().parent
rows=json.loads((R/'fresh_reference_results.json').read_text())
old=np.load(R/'reference_inputs/gb_eta010_g256_mean.npz')['G']
full=np.load(R/'shift_candidates.npz')['G'][-1]
quad=np.load(R/'candidate_frequency.npz')['G']
A=np.load(R/'fresh_G_N1536_eta0.005.npz')['G'];B=np.load(R/'fresh_G_N1536_eta0.0025.npz')['G'];G=2*B-A
IL=np.load(R/'fresh_residue_N1536.npz')['ImL'];I=(G-G.conj().T)/(2j);IT=I-IL
g=[.1,.2,.2,.2,.2,.1]
u=evaluate(1.2,g,G,{'qL':IL,'qT':IT});u0=evaluate(1.2,g,G,{'qL':IL,'qT':IT},eta_lead=0);v=p.solve_junction(g,1.2,G,IT,IL)
checks={'software':{'python':platform.python_version(),'numpy':np.__version__},'quadrant_vs_full_same_grid_max_abs':float(abs(full-quad).max()),'finite_core_same_inputs_two_implementations':{'qT_difference':float(abs(u0['partials']['qT']-v['TqT'])),'qL_difference':float(abs(u0['partials']['qL']-v['TqL'])),'reflection_difference':float(abs(u0['Rq']-v['R'])),'line_transmission_difference':float(abs(u0['Tq']-v['Tq'])),'eta_lead':0,'lead_eta_1e12_vs_zero_qT_abs':float(abs(u['partials']['qT']-u0['partials']['qT'])),'note':'Same kernels and physical model; a check of finite-core algebra, not independent bulk integration accuracy.'},'spectral_matrices_N1536':{'total_min_eig':float(np.linalg.eigvalsh(I).min()),'qL_min_eig':float(np.linalg.eigvalsh(IL).min()),'subtracted_qT_min_eig':float(np.linalg.eigvalsh(IT).min())},'low_pair_last_grid_change_abs':float(abs(rows[-1]['results']['low_pair_linear_richardson']['partials']['qT']-rows[-2]['results']['low_pair_linear_richardson']['partials']['qT'])),'N1536_fit_change_abs':float(abs(rows[-1]['results']['low_pair_linear_richardson']['partials']['qT']-rows[-1]['results']['three_point_quadratic']['partials']['qT'])),'fresh_minus_retained_central_qT':float(u['partials']['qT']-.1256063099667304),'limitations':['Regulator and shift here are chosen and explicit; this does not authenticate the missing historical producer or convention.','Grid/fit differences are sensitivities, not rigorous error bounds.','Fresh low-pair row is an additional corroboration; original table values remain unchanged.']}
(R/'summary_checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n')
print(json.dumps(checks,indent=2))
