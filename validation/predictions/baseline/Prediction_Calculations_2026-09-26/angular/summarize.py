"""Numerical summary and scientific plot; no additional physics calculation."""
from pathlib import Path
import csv,json,hashlib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
f=json.loads((P/'integral_fine.json').read_text());m=json.loads((P/'integral_medium.json').read_text())
s=json.loads((P/'integral_source_N160.json').read_text());angles=[.5,1,2,5,10,20,30]
rows=[]
for a in angles:
 c=np.cos(np.deg2rad(a));# stable positive formula for the point-force transverse integral
 ws=(1-c)**2*(c+2)/3;wl=(1-c)*(1+c+c*c)/(3*3**1.5)
 pointshare=ws/(ws+wl)
 for kind in ['group','wave']:
  for side in ['forward','backward']:
   key=f'{kind}_{side}_{a:g}'
   for name in ['full','X_only','R_only','same_X','same_R','same_cross']:
    x=f['cone_branch_powers'][key][name];q=sum(x[:2]);l=x[2]
    q0=sum(m['cone_branch_powers'][key][name][:2]);l0=m['cone_branch_powers'][key][name][2]
    row={'half_angle_deg':a,'direction_kind':kind,'side':side,'state':name,'incident_qT':q,'incident_qL':l,
         'qT_share':q/(q+l) if name!='same_cross' else None,
         'qT_refinement_relative':abs(q-q0)/max(abs(q),1e-300),'qL_refinement_relative':abs(l-l0)/max(abs(l),1e-300),
         'point_force_continuum_share':pointshare if side=='forward' and kind=='group' else None}
    if a in [10,20,30]:
     y=s['cone_branch_powers'][key][name]
     row['qT_N160_source_relative_change']=sum(y[:2])/q-1
     row['qL_N160_source_relative_change']=y[2]/l-1
    rows.append(row)
with (P/'cone_results.csv').open('w') as h:
 keys=sorted(set().union(*(r.keys() for r in rows)));w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)
overall={name:{'qT':sum(b[:2]),'qL':b[2],'bulk':sum(b),'share':sum(b[:2])/sum(b)} for name,b in f['total_branch_powers'].items() if name!='same_cross'}
summary={'status':'completed_numerical_prediction_with_explicit_limits','reference_omega':1.2,'global':overall,
         'rows':rows,'source_convergence':'N160 versus selected N192 imaginary source kernel; identical selected PV. Not a rigorous uncertainty enclosure.',
         'point_force_comparator':{'cT':1.,'cL':math.sqrt(3),'normalization':'Cone shear fraction only; prescribed force amplitude cancels. No incident channel or absolute incident-normalized power assigned.',
         'weights_up_to_common_2pi':{'S':'(1-cos(theta))^2*(2+cos(theta))/3','P':'(1-cos(theta)^3)/(3*(sqrt(3))^3)'}},
         'preimage_scope':'All radial roots globally complete by exact reference-symbol proof; all finite-cone polar intervals numerically searched on two azimuth grids. No certified exclusion of arbitrarily narrow unsampled angular islands; no exact far-field dark-cone claim.',
         'interpretation':'Ideal angular modal power collector in the asymptotic/incoherent channel sense. Not a coherent receiver, loaded detector, finite-distance field measurement, or new experimental validation.'}
(P/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')

plt.rcParams.update({'font.size':9,'axes.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(8.0,3.35),layout='constrained')
styles={'full':('Full','#0072B2','o','-'),'X_only':('X-only','#D55E00','s','--'),'R_only':('R-only','#009E73','^','-.')}
for name,(label,color,marker,line) in styles.items():
 d=[r for r in rows if r['direction_kind']=='group' and r['side']=='forward' and r['state']==name]
 ax[0].plot(angles,[100*r['qT_share'] for r in d],label=label,color=color,marker=marker,ls=line,ms=4)
 ax[1].plot(angles,[r['incident_qT'] for r in d],label=label,color=color,marker=marker,ls=line,ms=4)
point=[next(r['point_force_continuum_share'] for r in rows if r['half_angle_deg']==a and r['direction_kind']=='group' and r['side']=='forward') for a in angles]
ax[0].plot(angles,100*np.array(point),color='0.25',ls=':',marker='D',ms=3,label='Axial point force (continuum)')
ax[0].set_ylabel('qT / (qT + qL) in cone (%)');ax[0].set_ylim(-2,103)
ax[1].set_ylabel('Selected qT / incident power');ax[1].set_yscale('log')
for a in ax:
 a.set_xscale('log');a.set_xticks([.5,1,2,5,10,20,30],['0.5','1','2','5','10','20','30']);a.set_xlabel('Forward energy-flow cone half-angle (deg)');a.grid(alpha=.2)
ax[0].legend(fontsize=7,loc='center left',frameon=False)
ax[1].legend(fontsize=8,loc='upper left',frameon=False)
fig.suptitle(r'Reference junction, $\omega=1.2$: cone purity and collection strength',fontsize=10)
fig.savefig(P/'angular_collection.png',dpi=220,bbox_inches='tight')
plt.close(fig)
files={str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in P.rglob('*') if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in str(p)}
(P/'MANIFEST.json').write_text(json.dumps({'files':files,'original_archive_changed':False,'inputs':'Small matrix/metadata copies; no archive modules imported.'},indent=2)+'\n')
print(json.dumps({'global':overall,'output_rows':len(rows),'files':len(files)},indent=2))
