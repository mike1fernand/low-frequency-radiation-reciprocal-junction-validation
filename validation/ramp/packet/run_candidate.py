"""Bounded staged packet run; preserve all states/history and report measured cost."""
from pathlib import Path
import argparse, hashlib, json, os, resource, subprocess, time
import numpy as np

ROOT=Path(__file__).resolve().parent
ENV=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
G=[(j+1)/7 for j in range(6)]
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,x:p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def run(out,target,half=28,initial=None,resume=None,every=300):
    if not out.exists():out.mkdir()
    cmd=[str(ROOT/'evolve_candidate'),str(out),str(target),'.04',str(half)]+[format(x,'.17g') for x in G]+[str(every)]
    if resume:cmd.append(str(resume))
    begin=time.perf_counter()
    result=subprocess.run(cmd,input=initial,env=ENV,capture_output=True,text=True)
    (out/f'run_to_{target}.stdout').write_text(result.stdout+'\n'+result.stderr)
    if result.returncode:raise RuntimeError(result.stderr or result.stdout)
    d=json.loads(result.stdout);d['wall_seconds']=time.perf_counter()-begin
    d['max_rss_kib_process_tree']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    J(out/f'RUN_TO_{target}.json',d)
    return d

ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['pilot','full']);args=ap.parse_args()
if args.stage=='pilot':
    j=np.arange(120);test=np.r_[.003*np.sin(.17*j)+.001*np.cos(.31*j),.001*np.cos(.19*j)]
    test_text=' '.join(format(float(x),'.17g') for x in test)
    c=ROOT/'resume_control_continuous';s=ROOT/'resume_control_staged'
    assert not c.exists() and not s.exists()
    run(c,20,2,test_text,every=10);run(s,10,2,test_text,every=10);run(s,20,2,resume=s/'checkpoint_10.bin',every=10)
    for name in ['checkpoint_20.bin','line.bin','corner.bin','monitor.csv']:
        assert (c/name).read_bytes()==(s/name).read_bytes(),name
    J(ROOT/'RESUME_CONTROL.json',{'status':'PASS','dt':.04,'H':2,'staged_steps':[10,20],'continuous_steps':20,'comparison':'Fullcheckpoint and every line/corner/history byte match exactly; random nonzero line initialstate exercisesactualcontact.','source_sha256':H(ROOT/'evolve_candidate.cpp')})
    n=120;k=2*np.pi*np.fft.fftfreq(n);w=2.6*abs(np.sin(k/2));sense=np.sign(k);sense[60]=0
    j=np.arange(n)-18;q=.002*np.exp(-j*j/98)*np.cos(2*np.arcsin(6/13)*j)
    v=np.fft.ifft(-1j*sense*w*np.fft.fft(q)).real
    E0=float((v@v+1.69*np.diff(q)@np.diff(q))/2)
    Hc0=float(2/3*np.sum(np.array(G)*(q[60:66]**2+q[61:67]**2)))
    init=' '.join(format(float(x),'.17g') for x in np.r_[q,v]);(ROOT/'initial_qv.txt').write_text(init+'\n')
    np.savez_compressed(ROOT/'INITIAL_STATE.npz',q=q,v=v)
    meta={'N_line':120,'bulk_shape':[120,56,56,3],'quarter_shape':[120,28,28,3],'profile':G,'frequency_carrier':1.2,'carrier_k':float(2*np.arcsin(6/13)),'amplitude':.002,'gaussian_sigma':7.,'center_index':18,'contact_segments':[60,61,62,63,64,65],'dt':.04,'target_final_step':1425,'target_time':57.,'save_cadence':.2,'snapshot_start':48.,'initial_open_line_energy':E0,'initial_contact_energy':Hc0,'initial_contact_over_line':Hc0/E0,'scalar_line_stiffness':1.69,'masses':1.,'bulk_axis_and_face_stiffness':1.,'B_plus':[[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],'boundary':'Allphysicalouterboundariesfree;noaxialortransversewrapbonds. Onlyexact y/zreflectionreduction;noaxialreflections.','checkpoint':'Binarycompletequarterbulk+fullline x/v atintegralVerletstep, time,dt,profile,Bplus,constants,cumulativecurrent,lastJ,initialH; reconstructfullbulkexactlybyrecordedtransversereflections.','source_sha256':H(ROOT/'evolve_candidate.cpp'),'binary_sha256':H(ROOT/'evolve_candidate'),'initial_arrays_sha256':H(ROOT/'INITIAL_STATE.npz')}
    J(ROOT/'PREDECLARED_PACKET_CONFIGURATION.json',meta)
    out=ROOT/'candidate_dt004_box56';assert not out.exists()
    d=run(out,50,28,init)
    projected=d['assembly_seconds']+d['evolution_seconds']/50*1425
    J(ROOT/'PILOT_CHECK.json',{'status':'PASS','pilot':d,'projected_full_evolve_seconds':projected,'resource_cap_seconds':300,'resource_cap_rss_kib':2097152,'fullrun_within_measured_cost_gate':projected<=300 and d['max_rss_kib_process_tree']<=2097152,'no_stationary_target_used_in_selector':True})
    print(json.dumps({'pilot':d,'projected_full_evolve_seconds':projected,'initial_open_E0':E0,'initial_contact_E0':Hc0,'fullrun_gate':projected<=300 and d['max_rss_kib_process_tree']<=2097152},indent=2))
else:
    pilot=json.loads((ROOT/'PILOT_CHECK.json').read_text());assert pilot['fullrun_within_measured_cost_gate']
    out=ROOT/'candidate_dt004_box56';assert not (out/'checkpoint_1425.bin').exists()
    d=run(out,1425,resume=out/'checkpoint_50.bin')
    print(json.dumps(d,indent=2))
