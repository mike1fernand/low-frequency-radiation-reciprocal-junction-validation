"""Confirm screen candidates with retained refinements and independent tetra kernels."""
import csv,json
from pathlib import Path
import numpy as np
from interference_screen import HERE,ROOT,ARCHIVE,MAIN,sha,row,load_kernel,path_id,resolve_id,main_binding

profile=[.1,.2,.2,.2,.2,.1]
paths=[ARCHIVE/'surface/selected_nodes/N128/omega0.8.npz',
       ARCHIVE/'surface/selected_nodes/N160/omega0.8.npz',
       ARCHIVE/'tetra/packet_N96/omega_000.npz',
       ARCHIVE/'tetra/packet_N128/omega_000.npz']
reversals=[row(p,cq,profile) for p in paths for cq in [.8,1.3]]
reference=[row(ARCHIVE/'selected/REFERENCE_KERNEL.npz',1.3,profile),
           row(ARCHIVE/'tetra/reference_N192/omega_000.npz',1.3,profile)]

# Further alpha points only at confirmed candidate, on independent N160/N128 kernels.
from interference_screen import solve
alpha_refinement=[]
for p in [paths[1],paths[3]]:
    G,Q,m=load_kernel(p)
    for a in [.25,.75]:
        alpha_refinement.append({'kernel':path_id(p),'kernel_sha256':sha(p),
            'omega':.8,'cq':.8,'alpha':a,'result':solve(G,Q,.8,.8,profile,ar=a)})

# A complete motif means all four corners and both links to each corner, ns=1.
# The local 2x2 block of the translationally invariant bare kernel is exact for
# this support; it is not a truncation of the radiating exterior.
motif_paths=[paths[1],paths[3],ARCHIVE/'selected/REFERENCE_KERNEL.npz',
             ARCHIVE/'tetra/reference_N192/omega_000.npz',
             ARCHIVE/'surface/selected_nodes/N160/omega1.6.npz']
motifs=[]
for p in motif_paths:
    G,Q,m=load_kernel(p,1);P=np.diag([-1.,1.]);qt=Q[0]+Q[1]
    for cq in [1.,1.3,2.]:
        for g in [.1,.2,.4]:
            r=row(p,cq,[g])
            if r:
                r['reflection']={'P':P.tolist(),'G_commutator_max':float(abs(P@G-G@P).max()),
                    'QqT_commutator_max':float(abs(P@qt-qt@P).max()),'G_XR_abs':float(abs(G[0,1])),
                    'QqT_XR_abs':float(abs(qt[0,1]))}
                motifs.append(r)
two_motifs=[row(p,1.3,[.2,.2]) for p in [ARCHIVE/'selected/REFERENCE_KERNEL.npz',ARCHIVE/'tetra/reference_N192/omega_000.npz']]

allrows=reversals+reference+motifs+two_motifs
used={r['kernel'] for r in allrows}
bindings=[]
for rel in sorted(used):
    p=resolve_id(rel); meta=p.with_suffix('.json')
    bindings.append({'path':rel,'sha256':sha(p),'metadata_path':path_id(meta),'metadata_sha256':sha(meta)})
checks={
 'reversal_confirmed_all_four_kernels_both_speeds':all(r['redressed_full_minus_X']<0 for r in reversals),
 'same_state_cross_positive_all_reversal_rows':all(r['full']['same_state_qT']['XR']>0 for r in reversals),
 'reference_redressed_increase_both_methods':all(r['redressed_full_minus_X']>0 for r in reference),
 'single_motif_max_abs_XR':max(abs(r['full']['same_state_qT']['XR']) for r in motifs),
 'single_motif_max_relative_XR':max(abs(r['full']['same_state_qT']['normalized_cross']) for r in motifs),
 'single_motif_max_G_commutator':max(r['reflection']['G_commutator_max'] for r in motifs),
 'single_motif_max_Q_commutator':max(r['reflection']['QqT_commutator_max'] for r in motifs),
 'two_motif_cross_positive_both_methods':all(r['full']['same_state_qT']['XR']>0 for r in two_motifs),
 'max_full_flux_closure':max(r['full']['closure'] for r in allrows),
 'max_augmented_eta_difference':max(r['full']['augmented_eta_error'] for r in allrows),
 'max_paired_vs_X_plus_R_Hessian_difference':max(r['full']['hessian_pair_error'] for r in allrows),
 'max_same_state_decomposition_residual':max(r['full']['decomposition_residual'] for r in allrows),
}
out={'status':'confirmed_finite_numerical_and_analytic_results_not_interval_certified',
 'scope':'Independent Cartesian finite assembly; archived surface/PV and independent tetra bulk kernels, with no fresh bulk integration.',
 'main_source':main_binding(),
 'scripts':{p.name:sha(p) for p in [HERE/'interference_screen.py',Path(__file__)]},
 'inputs':bindings,'checks':checks,'reversal_rows':reversals,'reference_rows':reference,
 'alpha_refinement':alpha_refinement,'single_motif_rows':motifs,'two_motif_positive_controls':two_motifs}
(HERE/'confirmation.json').write_text(json.dumps(out,indent=2))
flat=[]
for kind,rows in [('reversal',reversals),('reference',reference),('single_motif',motifs),('two_motif',two_motifs)]:
    for r in rows:
        flat.append({'kind':kind,'omega':r['omega'],'cq':r['cq'],'profile':json.dumps(r['profile']),
          'method':'tetra' if '/tetra/' in r['kernel'] else 'surface_PV','grid':r['grid'],
          'TqT_alpha0':r['controls']['0.0']['TqT'],'TqT_alphaHalf':r['controls']['0.5']['TqT'],
          'TqT_alpha1':r['full']['TqT'],'same_state_XX':r['full']['same_state_qT']['XX'],
          'same_state_RR':r['full']['same_state_qT']['RR'],'same_state_XR':r['full']['same_state_qT']['XR'],
          'full_minus_X':r['redressed_full_minus_X'],'kernel_sha256':r['kernel_sha256']})
with (HERE/'confirmation.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
print(json.dumps(checks,indent=2))
