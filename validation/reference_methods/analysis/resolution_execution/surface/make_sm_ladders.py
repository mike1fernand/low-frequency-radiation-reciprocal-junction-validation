from pathlib import Path
import json,sys,hashlib,math
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(HERE/'src'));sys.path.insert(0,str(HERE))
import spectral_histogram as sh
from evaluate_refinement import transmit
paths=[HERE/'histograms'/f'sobol_seed{s}_p23_nseg10_bins16384.npz' for s in [20260920,20260921]]
hists=[sh.load_histogram(p) for p in paths];groups={M:[h.principal_value_jackson(1.44,M)[:12,:12] for h in hists] for M in [192,256,320]};means={M:np.mean(a,axis=0) for M,a in groups.items()};ref=means[256];pvrows=[]
for M in [192,256,320]:
 entries=[float(a[0,0]) for a in groups[M]];mean_entry=float(ref[0,0]) if M==256 else float(means[M][0,0]);pvrows.append({'M':M,'mean_PV_Frobenius_norm':float(np.linalg.norm(means[M])),'d_M_relative_to_M256':float(np.linalg.norm(means[M]-ref)/np.linalg.norm(ref)),'s00_mean':mean_entry,'s00_seed_min':min(entries),'s00_seed_max':max(entries),'s00_seed_span_over_mean_percent':100*(max(entries)-min(entries))/abs(mean_entry)})
shellrows=[]
for N in [128,160,192]:
 p=HERE/'results'/f'omega1.2_N{N}.npz';r=json.loads(p.with_suffix('.json').read_text());a=np.load(p);br=[a[f'Im{i}'] for i in range(3)];co=transmit(1.2,ref,br)
 shellrows.append({'N':N,'Im_total_Frobenius_norm':float(np.linalg.norm(a['Im_total'])),'Im_branch2_Frobenius_norm':float(np.linalg.norm(br[2])),'min_eigenvalue_Im_total_over_pi':float(np.linalg.eigvalsh(a['Im_total']/math.pi).min()),'raw_symmetry_residual':r['raw_reciprocity_residual'],'qT_same_selected_p23_M256':co['qT'],'matrix_path':str(p.relative_to(ROOT)),'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
out={'omega':1.2,'physical_block_shape':[12,12],'PV_definition':'equal mean of two separately scrambled Sobol sequences, seeds20260920 and20260921,2^23 samples per seed,16384 energybins on[0,12],margin1.02; table is first12×12 physicalblock','PV_rows':pvrows,'surface_rows':shellrows,'support_scope':'Histograms have10 sampled separationblocks; first12PV equivalence checked against prior6 support. Spectralnseg10 versus6 top12 equivalence checked onlyatN192.','producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'SM_MATRIX_LADDERS.json').write_text(json.dumps(out,indent=2)+'\n')
lines=['# Current physical12×12 matrix ladders','',out['PV_definition']+'.','', '| M | Mean PV Frobeniusnorm | Relative difference fromM256 | s00 seedspan/mean (%) |','|---:|---:|---:|---:|']
for r in pvrows:lines.append(f"|{r['M']}|{r['mean_PV_Frobenius_norm']:.9g}|{r['d_M_relative_to_M256']:.6g}|{r['s00_seed_span_over_mean_percent']:.6g}|")
lines+=['','| N | Total Im norm | Branch2 Im norm | min eig(Im/pi) | Raw symmetryresidual | qT usingselectedM256 |','|---:|---:|---:|---:|---:|---:|']
for r in shellrows:lines.append(f"|{r['N']}|{r['Im_total_Frobenius_norm']:.9g}|{r['Im_branch2_Frobenius_norm']:.9g}|{r['min_eigenvalue_Im_total_over_pi']:.6g}|{r['raw_symmetry_residual']:.6g}|{r['qT_same_selected_p23_M256']:.10f}|")
lines+=['',out['support_scope'],'','Observed sensitivity diagnostics; no exact-lattice errorbound.']
(HERE/'SM_MATRIX_LADDERS.md').write_text('\n'.join(lines)+'\n');print(json.dumps(out,indent=2))
