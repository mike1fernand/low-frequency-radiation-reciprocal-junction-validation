#!/usr/bin/env python3
"""Fresh checkpointed spectral histograms with recorded deterministic prefixes.
The same scrambled Sobol prefix feeds increasing counts; these are dependent
refinements, while seeds are distinct. No statistical confidence bound is claimed.
"""
from pathlib import Path
import sys,json,time,argparse,hashlib,resource
import numpy as np
from scipy.stats import qmc
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'src'))
import spectral_histogram as sh
from model_api_v1_1 import ModelParameters
ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--max-power',type=int,default=22);ap.add_argument('--bins',type=int,default=16384);ap.add_argument('--nseg',type=int,default=6);ap.add_argument('--resume-power',type=int);a=ap.parse_args()
if a.nseg not in [6,10]:raise ValueError('Only declared6/10 support lengths are admitted')
if not 20<=a.max_power<=23:raise ValueError('Bounded worker supports checkpoint powers20–23 only')
if a.bins not in [8192,16384]:raise ValueError('Unsupported spectral bin count')
p=ModelParameters();edges=np.linspace(0.,12.,a.bins+1);rawmass=np.zeros((a.nseg,2,2,a.bins),complex);engine=qmc.Sobol(3,scramble=True,seed=a.seed)
start=time.time();stages=[];max_gap=0.;tiecount=0;split_bin_count=0;chunk=65536;sample_hash=hashlib.sha256()
resume_count=0;resume_source=None
if a.resume_power is not None:
 if a.resume_power!=22 or a.max_power!=23 or a.nseg!=10:raise ValueError('Only validatedp22→p23 nseg10 extension allowed')
 resume_source=HERE/'histograms'/f'sobol_seed{a.seed}_p22_nseg10_bins{a.bins}.npz'
 prior=sh.load_histogram(resume_source)
 prior_record=json.loads(resume_source.with_suffix('.json').read_text())
 if hashlib.sha256(resume_source.read_bytes()).hexdigest()!=prior_record['histogram_sha256']:raise RuntimeError('Resume histogramhash mismatch')
 if prior.sample_count!=2**22 or prior.nseg!=a.nseg or prior.metadata['seed']!=a.seed or not np.array_equal(prior.edges,edges):raise RuntimeError('Resume configuration mismatch')
 resume_count=prior.sample_count;rawmass[:]=prior.mass*resume_count*p.m_bulk
# _acc's per-point weights are accumulated at norm1; checkpoints normalize by m*n.
for st in range(0,2**a.max_power,chunk):
 points=engine.random(chunk);sample_hash.update(points.tobytes())
 if st+chunk<=resume_count:
  if st+chunk==resume_count and sample_hash.hexdigest()!=prior.metadata['sample_sequence_sha256']:raise RuntimeError('Resumed Sobolprefix mismatch')
  continue
 k=2*np.pi*(points-.5)
 vals,c=sh._mode_data(k[:,0],k[:,1],k[:,2],p)
 if not np.all(np.isfinite(vals)) or not np.all(np.isfinite(c)):raise RuntimeError('Nonfinite spectral sample')
 if vals.min() < -256*np.finfo(float).eps or vals.max()>12+256*np.finfo(float).eps:raise RuntimeError('Spectrum exceeds declared bin range')
 tol=64*np.finfo(float).eps*np.maximum(1.,np.max(abs(vals),axis=-1))
 ties=np.diff(vals,axis=-1)<=tol[:,None];tiecount+=int(ties.sum())
 idx=np.clip(np.searchsorted(edges,vals,side='right')-1,0,a.bins-1)
 split_bin_count+=int((ties & (np.diff(idx,axis=-1)!=0)).sum())
 if split_bin_count:raise RuntimeError('Floating-degenerate spectral pair split across histogram bins; explicit treatment required')
 sh._acc(rawmass,vals,c,k[:,0],edges,a.nseg,1.)
 count=st+chunk
 if count not in [2**x for x in range(20,a.max_power+1)]:continue
 power=count.bit_length()-1
 metadata={'seed':a.seed,'power':power,'nbins':a.bins,'energy_range':[0.,12.],'nseg':a.nseg,'chunk':chunk,'sequence':'SciPy scrambled Sobol random(65536), contiguous initial prefix','sample_sequence_sha256':sample_hash.hexdigest(),'prefixes_are_nested_not_independent':True,'resumed_from_count':resume_count,'resumed_from_histogram_sha256':None if resume_source is None else hashlib.sha256(resume_source.read_bytes()).hexdigest()}
 h=sh.SpectralHistogram(edges.copy(),rawmass.copy()/(p.m_bulk*count),count,'sobol',metadata)
 support='' if a.nseg==6 else f'_nseg{a.nseg}'
 tag=f'sobol_seed{a.seed}_p{power}{support}_bins{a.bins}';file=HERE/'histograms'/(tag+'.npz');sh.save_histogram(file,h)
 expected=np.zeros((a.nseg,2,2),complex);expected[0]=np.eye(2)
 residual=float(abs(h.mass.sum(axis=-1)-expected).max())
 if residual>1e-4:raise RuntimeError('Gross sample mass sum-rule failure')
 rec={'schema':'paper1.fresh-spectral-histogram.v1','sampler':'sobol','sample_count':count,'elapsed_seconds':time.time()-start,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'metadata':metadata,'sum_rule_max_abs_error':residual,'degenerate_pair_samples':tiecount,'degenerate_pairs_split_across_bins':split_bin_count,'histogram_file':file.name,'histogram_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'spectral_source_sha256':hashlib.sha256((HERE/'src/spectral_histogram.py').read_bytes()).hexdigest()}
 (HERE/'histograms'/(tag+'.json')).write_text(json.dumps(rec,indent=2)+'\n');stages.append(rec);print(json.dumps(rec),flush=True)
(HERE/'histograms'/(f'seed{a.seed}'+('' if a.nseg==6 else f'_nseg{a.nseg}')+('_STAGES.json' if not resume_count else '_p23_EXTENSION.json'))).write_text(json.dumps(stages,indent=2)+'\n')
