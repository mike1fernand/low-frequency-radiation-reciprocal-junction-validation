"""Project a freshly generated reference trajectory and compare with archived observables."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import struct

import numpy as np

from verify_packet_reference import DATA, check, incident_measure, mean, read

# Reference crystal in a=m_b=K_axis=K_diag=1 units; one direction per undirected bond.
BONDS = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),
         (1,0,-1),(0,1,1),(0,1,-1)]


def eigenmodes(shape):
    k = np.meshgrid(*[2*np.pi*np.fft.fftfreq(n) for n in shape],indexing='ij',sparse=True)
    dynamical = np.zeros((*shape,3,3))
    for bond in BONDS:
        phase = sum(kj*bj for kj,bj in zip(k,bond))
        v = np.array(bond)
        # Sum the positive/negative neighbors' central longitudinal spring terms.
        dynamical += (4*np.sin(phase/2)**2)[...,None,None] * np.outer(v,v)/np.dot(v,v)
    return np.linalg.eigh(dynamical)


def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source,'sha256').hexdigest()


def read_field(stream, shape):
    n = int(np.prod(shape))*3
    raw = stream.read(4*n)
    if len(raw)!=4*n:
        raise ValueError("Incomplete field in snapshot")
    field = np.frombuffer(raw,dtype='<f4').reshape(*shape,3).astype('float64')
    return np.fft.fftn(field,axes=(0,1,2))


def analyze(path, trajectory_csv, checkpoint):
    with trajectory_csv.open() as current_csv, (DATA/'REFERENCE_TRAJECTORY_ROWS.csv').open() as stored_csv:
        new, old = list(csv.DictReader(current_csv)),list(csv.DictReader(stored_csv))
    fields=list(old[0])
    assert len(new)==len(old)==46 and list(new[0])==fields
    for key in fields:
        # The cutoff for the data files is checked separately from cross-build rounding differences.
        check([float(r[key]) for r in new],[float(r[key]) for r in old],2e-10)
    archive=read('REFERENCE_MODAL_ANALYSIS.json')
    _,_,e0,_,_,_=incident_measure()
    reference=archive['rows']
    selected=[]; max_periodic_difference=0.0; max_bound_difference=0.0
    with path.open('rb') as stream:
        count,*shape=struct.unpack('<QQQQ',stream.read(32))
        assert count==len(reference) and shape==[120,56,56]
        values,vectors=eigenmodes(shape)
        gap_tol=1.e-12*np.maximum(1.,np.max(abs(values),axis=-1))
        upper_ties=(values[...,2]-values[...,1])<=gap_tol
        is_triple=upper_ties & ((values[...,1]-values[...,0])<=gap_tol)
        double=upper_ties & ~is_triple
        assert (int(double.sum()),int(is_triple.sum()))==(132,231)
        norm=2*np.prod(shape)*e0
        for original in reference:
            t,=struct.unpack('<d',stream.read(8));check(t,original['time'],1e-11)
            u,v=read_field(stream,shape),read_field(stream,shape)
            cu=np.einsum('...ij,...i->...j',vectors,u)
            cv=np.einsum('...ij,...i->...j',vectors,v)
            energies=(abs(cv)**2+np.maximum(values,0)*abs(cu)**2)/norm
            periodic=float(np.sum(energies))
            max_periodic_difference=max(max_periodic_difference,abs(periodic-original['periodic_total_over_E0']))
            if 52.8-1e-8<=t<=55.4+1e-8:
                fixed=energies[...,2][~upper_ties].sum()
                ql_lower=ql_upper=fixed
                for mask,start in [(double,1),(is_triple,0)]:
                    v_block=cv[mask,start:];u_block=cu[mask,start:]
                    energy=(v_block[:,:,None]*v_block[:,None,:].conj()
                            +values[mask,start,None,None]*u_block[:,:,None]*u_block[:,None,:].conj())/norm
                    spectrum=np.linalg.eigvalsh(energy)
                    ql_lower+=spectrum[:,0].sum()
                    ql_upper+=spectrum[:,-1].sum()
                low,high=float(periodic-ql_upper),float(periodic-ql_lower)
                max_bound_difference=max(max_bound_difference,abs(low-original['qT_min_over_E0']),
                                         abs(high-original['qT_max_over_E0']))
                selected.append((t,low,high))
        assert stream.read(1)==b''
    t,lo,hi=np.array(selected).T
    low,high=mean(t,lo),mean(t,hi)
    assert len(selected)==14 and max_periodic_difference<3e-13 and max_bound_difference<3e-13
    digests={name:sha(p) for name,p in [('rows',trajectory_csv),('snapshots',path),('final_checkpoint',checkpoint)]}
    published={'rows':'e045b76ce94b5b033954243283638d3c191a4fb0f0dfe3bf8f61826221309cec',
               'snapshots':'920a0849ac1c787feb878ecc6c063b7da28b1f1343da0670af3c80bcc46e34da',
               'final_checkpoint':'b7da5019c405fee76dc7b56ed1c3a7056b5e1bf31b1356d39af105c5b2a2fa42'}
    return {'status':'PASS','integration_steps':1425,'dimensions':shape,'stored_snapshots':count,
            'reporting_samples':len(selected),'fresh_integrator':'archived single-thread velocity-Verlet implementation',
            'fresh_projection':'Explicit nine-family stiffness matrix; full symmetric eigensolve at every wavevector, degenerate Hermitian energy extrema',
            'qT_window_allocation_bounds':[low,high],'qT_midpoint':(low+high)/2,
            'max_periodic_bulk_deviation_from_archived_per_snapshot':max_periodic_difference,
            'max_qT_bound_deviation_from_archived_per_reporting_snapshot':max_bound_difference,
            'raw_file_sha256':digests,'all_raw_hashes_match_original_release':digests==published,
            'original_recipe':'PAPER1_VALIDATION_INTEGRATED_R8.zip / legacy_rc4/validation/direct_dynamics/run/PASS12_PRODUCTION_REGENERATION_RUNBOOK.md'}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('snapshots',type=Path);ap.add_argument('rows',type=Path);ap.add_argument('checkpoint',type=Path)
    args=ap.parse_args();print(json.dumps(analyze(args.snapshots,args.rows,args.checkpoint),indent=2))
