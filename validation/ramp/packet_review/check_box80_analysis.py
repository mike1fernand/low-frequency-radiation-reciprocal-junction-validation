"""Independent selector, normalization and three-snapshot modal checks.

Uses the exact closed bulk-symbol formula and total-minus-upper allocation,
rather than importing the packet analyzer. No new dynamics is performed.
"""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent;PACKET=ROOT.parent/'packet'
rp=PACKET/'CANDIDATE_PACKET_BOX80_RESULT.json';result=json.loads(rp.read_text());cfg=result['config']
data=PACKET/'candidate_dt004_box80';HSH=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=result['rows'];required=[]
for r in rows:
    passes=(r['link']<=.0025 and abs(r['J'])<=.0005 and r['shell']<=.002
            and abs(r['periodic_difference'])<=.0025 and r['drift']<=1e-5)
    assert passes==r['passed']
    if passes:required.append(r['t'])
runs=[]
for t in required:
    if not runs or t-runs[-1][-1]>.20000000001:runs.append([])
    runs[-1].append(t)
valid=[r for r in runs if r[-1]-r[0]>=2-1e-12]
best=max(valid,key=lambda x:(x[-1]-x[0],-x[0]))
assert [best[0],best[-1]]==result['longest']
assert result['selection_mode']=='coefficient_blind_search'

N=cfg['N_line'];M=cfg['bulk_shape'][1];h=M//2;n3=N*M*M
times=np.fromfile(data/'snap_times.bin',dtype=np.float64)
raw=np.memmap(data/'snap.bin',dtype=np.float64,mode='r').reshape(len(times),2,N+3*N*h*h)
line=np.memmap(data/'line.bin',dtype=np.float64,mode='r').reshape(-1,2,N)
q0,v0=line[0];E0=.5*(v0@v0+cfg['scalar_line_stiffness']*(np.diff(q0)@np.diff(q0)))
assert abs(E0-result['open_line_E0'])<1e-20
assert np.allclose(cfg['profile'],np.arange(1,7)/7,rtol=0,atol=0)
assert cfg['contact_segments']==list(range(60,66))

kx=(2*np.pi*np.fft.fftfreq(N))[:,None,None]
ky=(2*np.pi*np.fft.fftfreq(M))[None,:,None]
kz=(2*np.pi*np.fft.rfftfreq(M))[None,None,:]
cx,cy,cz=np.cos(kx),np.cos(ky),np.cos(kz)
sx,sy,sz=np.sin(kx),np.sin(ky),np.sin(kz)
D=np.zeros((N,M,h+1,3,3))
D[...,0,0]=2*(3-cx*(1+cy+cz));D[...,1,1]=2*(3-cy*(1+cx+cz));D[...,2,2]=2*(3-cz*(1+cx+cy))
D[...,0,1]=D[...,1,0]=2*sx*sy;D[...,0,2]=D[...,2,0]=2*sx*sz;D[...,1,2]=D[...,2,1]=2*sy*sz
lam,vec=np.linalg.eigh(D);lam=np.maximum(lam,0)
w=np.broadcast_to(np.r_[1.,np.full(h-1,2.),1.],lam.shape[:-1])
top=(lam[...,2]-lam[...,1])<=1e-12*np.maximum(1,lam[...,2])
triple=top & ((lam[...,2]-lam[...,0])<=1e-12*np.maximum(1,lam[...,2]))
assert int(np.sum(w[top & ~triple]))==result['upper_middle_doublet_points']
assert int(np.sum(w[triple]))==result['triple_points']
def unfold(row):
    q=row[N:].reshape(N,h,h,3)
    z=np.concatenate([q[:,:,::-1]*[1,1,-1],q],axis=2)
    return np.concatenate([z[:,::-1]*[1,-1,1],z],axis=1)
samples=np.asarray(result['branch_samples'])
chosen=[0,len(samples)//2,len(samples)-1];checks=[]
for j in chosen:
    t=samples[j,0];index=np.argmin(abs(times-t));assert abs(times[index]-t)<1e-12
    x=unfold(raw[index,0]);v=unfold(raw[index,1])
    uh=np.fft.rfftn(x,axes=(0,1,2));vh=np.fft.rfftn(v,axes=(0,1,2))
    kinetic=np.sum(abs(vh)**2,axis=-1)
    potential=np.einsum('...i,...ij,...j->...',uh.conj(),D,uh).real
    factor=w/(2*n3*E0);total=float(np.sum(factor*(kinetic+potential)))
    # Upper-separated qT is total minus the one-dimensional upper energy.
    upper=np.einsum('...i,...i->...',vec[...,2],uh)
    upper_v=np.einsum('...i,...i->...',vec[...,2],vh)
    fixed=float(np.sum((factor*(kinetic+potential-(abs(upper_v)**2+lam[...,2]*abs(upper)**2)))[~top]))
    # At top ties, allocate the full energy-coherence matrix invariantly.
    pu=np.einsum('...ia,...i->...a',vec,uh);pv=np.einsum('...ia,...i->...a',vec,vh)
    first_energy=factor*(abs(pv[...,0])**2+lam[...,0]*abs(pu[...,0])**2)
    fixed+=float(np.sum(first_energy[top & ~triple]));lo=hi=fixed
    for start,rank,mask in [(1,1,top & ~triple),(0,2,triple)]:
        a=pu[mask,start:];b=pv[mask,start:];omega2=lam[mask,start:].mean(axis=-1)
        gram=(b[:,:,None]*b[:,None,:].conj()+omega2[:,None,None]*a[:,:,None]*a[:,None,:].conj())*factor[mask,None,None]
        eigen=np.linalg.eigvalsh(gram)
        lo+=float(eigen[:,:rank].sum());hi+=float(eigen[:,-rank:].sum())
    got=np.array([lo,hi,total]);error=float(np.max(abs(got-samples[j,1:])))
    assert error<1e-12,(t,error,got,samples[j])
    parseval=float(abs(np.sum(abs(vh)**2*w[...,None])/n3-np.sum(v*v))/max(np.sum(v*v),1e-30))
    assert parseval<1e-12
    checks.append({'time':float(t),'independent_low_high_total':got.tolist(),'stored_low_high_total':samples[j,1:].tolist(),
                   'maximum_abs_difference':error,'kinetic_parseval_relative_difference':parseval})

mean=np.trapezoid(samples[:,1:],x=samples[:,0],axis=0)/(samples[-1,0]-samples[0,0])
assert np.max(abs(mean-np.asarray(result['mean_low_high_total'])))<1e-15
old_coarse=samples[::3];old_window=[float(old_coarse[0,0]),float(old_coarse[-1,0])]
indices=np.unique(np.r_[np.arange(0,len(samples),3),len(samples)-1]);coarse=samples[indices]
same_window_mean=np.trapezoid(coarse[:,1:],x=coarse[:,0],axis=0)/(coarse[-1,0]-coarse[0,0])
out={'status':'PASS','source_result_sha256_at_review':HSH(rp),'script_sha256':HSH(Path(__file__)),
     'independent_selector_window':[best[0],best[-1]],'all_five_predicates_recomputed':True,
     'selected_cadence_samples':len(samples),'three_snapshot_checks':checks,
     'primary_mean_recomputed':mean.tolist(),'legacy_stride3_window':old_window,
     'same_endpoint_stride3_window':[float(coarse[0,0]),float(coarse[-1,0])],
     'same_endpoint_coarse_mean':same_window_mean.tolist(),
     'cadence_erratum':'Original B[::3] omitted the final endpoint. It must not be labeled a same-window cadence comparison. Primary result unaffected; correction needs no new dynamics.',
     'scope':'Independent physical/FFT/selector checks of existing finite-packet results, not new long dynamics or an asymptotic branch-flux proof.'}
(ROOT/'BOX80_ANALYSIS_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
