"""Independent full-physical-spring one-step controls of reduced packet engine.

Two tiny transverse boxes; no long packet simulation. Complete physical forces
are assembled from undirected bonds, then compared to the actual binary's
velocity-Verlet step through a valid restart checkpoint.
"""
from pathlib import Path
import hashlib,json,os,struct,subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent
PACKET=ROOT.parent/'packet';N=120;dt=.04
g=np.arange(1,7,dtype=float)/7
B=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],dtype=np.int32)
HSH=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def unfold(v,H):
    a=v[N:].reshape(N,H,H,3)
    z=np.concatenate([a[:,:,::-1]*[1,1,-1],a],axis=2)
    return np.concatenate([z[:,::-1]*[1,-1,1],z],axis=1)
def physics(x,v,H):
    U=unfold(x,H);V=unfold(v,H);F=np.zeros_like(U);fq=np.zeros(N);pot=0.
    dq=np.diff(x[:N]);fq[:-1]+=1.69*dq;fq[1:]-=1.69*dq;pot+=.5*1.69*(dq@dq)
    shape=U.shape[:3]
    for b in B:
        low=tuple(slice(max(-d,0),s-max(d,0)) for d,s in zip(b,shape))
        high=tuple(slice(max(d,0),s-max(-d,0)) for d,s in zip(b,shape))
        delta=np.einsum('...c,c->...',U[high]-U[low],b)
        force=delta[...,None]*b/(b@b)
        F[low]+=force;F[high]-=force;pot+=float(np.sum(delta*delta)/(2*(b@b)))
    J=0.
    for j in range(6):
        i=60+j
        for e in [0,1]:
            sx=1-2*e
            for sy in [-1,1]:
                for sz in [-1,1]:
                    y=H-1 if sy<0 else H;z=H-1 if sz<0 else H
                    b=np.array([sx,sy,sz]);ext=b@U[i,y,z]-sx*x[i+e]
                    force=-g[j]*ext*b/3
                    F[i,y,z]+=force;fq[i+e]-=force[0]
                    pot+=g[j]*ext*ext/6
                    # Force on the line, paired with the sum of endpoint velocities.
                    J-=(-force)@(V[i,y,z]+np.array([v[i+e],0.,0.]))/2
    acceleration=np.r_[fq,F[:,H:,H:,:].reshape(-1)]
    reflected_force=unfold(acceleration,H)
    assert np.max(abs(reflected_force-F))<1e-15
    kinetic=.5*(v[:N]@v[:N]+np.sum(V*V))
    return acceleration,float(pot+kinetic),float(J)
def write_checkpoint(path,H,x,v,E,J):
    L=len(x)
    with path.open('wb') as f:
        f.write(struct.pack('<8sqiiid',b'P1PKT01\0',0,N,H,L,dt))
        f.write(g.astype('<f8').tobytes());f.write(B.astype('<i4').tobytes())
        f.write(np.array([1.,1.,1.,1.69],'<f8').tobytes())
        f.write(struct.pack('<ddd',0.,J,E));f.write(x.astype('<f8').tobytes());f.write(v.astype('<f8').tobytes())
def read_checkpoint(path,H):
    data=path.read_bytes();head=struct.unpack_from('<8sqiiid',data)
    assert head[:5]==(b'P1PKT01\0',1,N,H,N+3*N*H*H)
    assert head[5]==dt
    off=36+48+108+32
    cumulative,J,E=struct.unpack_from('<ddd',data,off);off+=24
    L=head[4];x=np.frombuffer(data,dtype='<f8',count=L,offset=off).copy();v=np.frombuffer(data,dtype='<f8',count=L,offset=off+8*L).copy()
    assert off+16*L==len(data)
    return x,v,cumulative,J,E
records=[]
for H,seed in [(2,20261108),(3,20261109)]:
    out=ROOT/f'random_step_H{H}';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(seed);L=N+3*N*H*H
    x=rng.normal(size=L)*.001;v=rng.normal(size=L)*.001
    a0,E0,J0=physics(x,v,H);vh=v+.5*dt*a0;expected_x=x+dt*vh
    a1,E1_unused,J1_unused=physics(expected_x,vh,H);expected_v=vh+.5*dt*a1
    _,E1,J1=physics(expected_x,expected_v,H)
    checkpoint=out/'independent_initial.bin';write_checkpoint(checkpoint,H,x,v,E0,J0)
    cmd=[str(PACKET/'evolve_candidate'),str(out),'1',str(dt),str(H),*[format(t,'.17g') for t in g],'1',str(checkpoint)]
    run=subprocess.run(cmd,capture_output=True,text=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
    assert run.returncode==0,run.stderr+run.stdout
    (out/'engine.stdout').write_text(run.stdout+run.stderr)
    gotx,gotv,cumulative,gotJ,gotE0=read_checkpoint(out/'checkpoint_1.bin',H)
    errors={'x':float(np.max(abs(gotx-expected_x))),'v':float(np.max(abs(gotv-expected_v))),
            'current':abs(gotJ-J1),'integrated_current':abs(cumulative-dt*(J0+J1)/2),'initial_energy':abs(gotE0-E0)}
    assert max(errors.values())<2e-15,errors
    records.append({'H':H,'seed':seed,'dofs':L,'errors':errors,'engine':json.loads(run.stdout),
                    'independent_physical_energy_initial':E0,'independent_physical_energy_final':E1})
result={'status':'PASS','scope':'Two one-step arbitrary-state controls of physical force/velocity-Verlet/current and y/z-only reflection reduction; not another long packet run.',
        'source_sha256':HSH(PACKET/'evolve_candidate.cpp'),'binary_sha256':HSH(PACKET/'evolve_candidate'),
        'script_sha256':HSH(Path(__file__)),'records':records}
(ROOT/'PHYSICAL_STEP_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
