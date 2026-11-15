#!/usr/bin/env python3
"""Checkpointed rigorous qT/Green integration on any one full angular tile.

Each accepted leaf contains its complete real parameter box. Pending leaves
are explicitly retained; a partial checkpoint is never a tile enclosure.
"""
import argparse, hashlib, json, math, os, sys, time
from fractions import Fraction
from pathlib import Path
from flint import arb,acb,ctx
from certified_tensor_sharp import gauss_box,NeighborhoodFailure
from source_moments import moment_values,symbol_matrix,coordinate_values,ball_serialization

def exact(s):
    q=Fraction(s);return arb(q.numerator)/q.denominator

def decode_ball(value):
    if 'real' in value:return acb(decode_ball(value['real']),decode_ball(value['imag']))
    def end(k):
        x=value[k];return arb(int(x['mantissa']))*arb(2)**int(x['exponent2'])
    return end('lower').union(end('upper'))

def source_hashes():
    return {name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('arbitrary_tile_sharp.py','certified_tensor_sharp.py','source_moments.py')}

def endpoint_certificate(bounds,chart,candidates=('1/2','9/16','5/8','11/16')):
    """Use a tile-wide Rayleigh bound for the top band; else global separator.

    Candidates are exact rational radii strictly below s=acos(3/4).
    lambda_2<=2(2-cos(t)-cos(t)^2), by rank-one interlacing.
    """
    a,b,c,d=map(exact,bounds);u=a.union(b);v=c.union(d)
    if not (0<=a and a<b and b<=1 and 0<=c and c<d and d<=1):raise ValueError('angular tile outside unit square')
    if chart not in (0,1,2):raise ValueError('invalid chart')
    s=(arb(3)/4).acos()
    for text in candidates:
        t=exact(text)
        if not (0<t and t<s):raise ValueError('candidate radius outside proven interval')
        D=symbol_matrix(t,u,v,chart)
        w=coordinate_values(arb(1),(a+b)/2,(c+d)/2,chart)
        ray=sum(w[i]*D[i,j]*w[j] for i in range(3) for j in range(3))/sum(x*x for x in w)
        lo=arb(36)/25-2*(2-t.cos()-t.cos()**2);hi=ray.real-arb(36)/25
        if lo>0 and hi>0 and ray.imag.contains(0):
            return t,{'kind':'tile_rayleigh','radius_rational':text,'bounds_rational':list(bounds),'chart':chart,
                      'middle_margin':ball_serialization(lo),'top_margin':ball_serialization(hi)}
    return s,{'kind':'global_separator','radius':'acos(3/4)','middle_margin':'13/200','top_margin':'3/50'}

def contour_legs(mode,height,t0=None):
    h=exact(height);R=(arb(1)/4).acos()
    if not h>0:raise ValueError('positive contour height required')
    if mode=='qt':
        if t0 is None:t0=(arb(3)/4).acos()
        return [('down',acb(t0),acb(0,-h)),('across',acb(t0,-h),acb(R-t0)),('up',acb(R,-h),acb(0,h))]
    if mode=='green':
        return [('down',acb(0),acb(0,-h)),('across',acb(0,-h),acb(R)),('up',acb(R,-h),acb(0,h)),('exterior',acb(R),acb(arb.pi()-R))]
    raise ValueError('mode must be qt or green')

def split_node(node,axis):
    lo,hi=map(Fraction,node['bounds'][axis]);mid=(lo+hi)/2
    children=[]
    for bit,ends in (('0',(lo,mid)),('1',(mid,hi))):
        b=[list(x) for x in node['bounds']];b[axis]=list(map(str,ends))
        children.append({'id':node['id']+bit,'leg':node['leg'],'bounds':b,'depth':node['depth']+1})
    return children

def audit_partition(state):
    """Rebuild every split exactly; require each generated leaf exactly once."""
    settings=state['settings'];mode=settings['mode']
    if mode not in ('qt','green'):raise ValueError('unknown contour mode')
    expected_legs=3 if mode=='qt' else 4
    if state['leg_count']!=expected_legs:raise ValueError('contour leg count inconsistent with mode')
    a,b,c,d=settings['bounds'];roots={f'{j}:':{'id':f'{j}:','leg':j,'bounds':[['0','1'],[a,b],[c,d]],'depth':0} for j in range(expected_legs)}
    expected=dict(roots)
    for split in state['splits']:
        node=expected.pop(split['id']);axis=split['axis']
        if axis not in (0,1,2):raise ValueError('invalid split axis')
        expected.update({x['id']:x for x in split_node(node,axis)})
    leaves=state['pending']+[v['node'] for v in state['accepted']]
    if len(leaves)!=len(expected):raise ValueError('coverage leaf count mismatch')
    observed={}
    for leaf in leaves:
        if leaf['id'] in observed:raise ValueError('duplicate leaf')
        if leaf!=expected.get(leaf['id']):raise ValueError('leaf geometry mismatch')
        observed[leaf['id']]=leaf
    if observed.keys()!=expected.keys():raise ValueError('missing coverage leaves')
    return len(leaves)

def atomic_json(dest,state):
    dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_suffix(dest.suffix+'.partial');tmp.write_text(json.dumps(state,indent=2)+'\n');os.replace(tmp,dest)

def recompute_total(state):
    total=[acb(0) for _ in range(17)]
    for leaf in state['accepted']:
        if len(leaf['moments'])!=17:raise ValueError('wrong moment count')
        for m,v in enumerate(leaf['moments']):
            value=decode_ball(v)
            if not value.is_finite():raise ValueError('nonfinite saved leaf')
            total[m]+=value
    return total

def completed_moments(state):
    if state['pending']:raise ValueError('complete claim has pending boxes')
    total=recompute_total(state)
    return [x.imag for x in total] if state['settings']['mode']=='qt' else total

def integrate_stage(settings,destination,max_seconds=20):
    if not math.isfinite(max_seconds) or max_seconds<=0:raise ValueError('finite positive time budget required')
    ctx.prec=settings['bits'];settings=dict(settings);dest=Path(destination)
    for key in ('bounds','orders'):settings[key]=list(settings[key])
    if len(settings['bounds'])!=4:raise ValueError('four angular bounds required')
    settings['bounds']=[str(Fraction(x)) for x in settings['bounds']]
    t0,cert=endpoint_certificate(settings['bounds'],settings['chart'],settings.get('endpoint_candidates',('1/2','9/16','5/8','11/16')))
    legs=contour_legs(settings['mode'],settings['height'],t0)
    identity={'sources':source_hashes(),'python':sys.version,'flint_version':__import__('flint').__version__}
    if dest.exists():
        state=json.loads(dest.read_text())
        if state['settings']!=settings or state['identity']!=identity:raise ValueError('checkpoint source/runtime/settings mismatch')
        if state['endpoint_certificate']!=cert:raise ValueError('endpoint certificate differs from recomputation')
        audit_partition(state)
        if state['status']=='rigorous_complete_tile_enclosure':
            if state['moments']!=list(map(ball_serialization,completed_moments(state))):raise ValueError('completed aggregate disagrees with saved leaves')
            return state
    else:
        a,b,c,d=settings['bounds']
        state={'status':'partial_checkpoint_no_tile_enclosure','settings':settings,'identity':identity,'endpoint_certificate':cert,'leg_count':len(legs),
               'pending':[{'id':f'{j}:','leg':j,'bounds':[['0','1'],[a,b],[c,d]],'depth':0} for j in range(len(legs))],
               'accepted':[],'splits':[],'stages':[]}
    start=time.monotonic();calls=0;accepted_before=len(state['accepted']);splits_before=len(state['splits'])
    state['status']='partial_checkpoint_no_tile_enclosure';state.pop('blocked_reason',None)
    density=exact(settings['density_tolerance'])/len(legs)
    while state['pending'] and time.monotonic()-start<max_seconds:
        node=state['pending'].pop(0);_,origin,delta=legs[node['leg']]
        def callback(q,u,v):
            nonlocal calls
            if time.monotonic()-start>=max_seconds:raise TimeoutError('stage time limit')
            calls+=1
            return [x*delta for x in moment_values(origin+delta*q,u,v,settings['chart'])]
        try:
            vals,meta=gauss_box(callback,[(exact(lo),exact(hi)) for lo,hi in node['bounds']],tuple(settings['orders']),tuple(settings.get('radius_multipliers',(2,4,8,16))),density)
            state['accepted'].append({'node':node,'moments':list(map(ball_serialization,vals)),'metadata':meta})
        except TimeoutError:
            state['pending'].insert(0,node);break
        except NeighborhoodFailure as exc:
            if node['depth']>=settings.get('max_depth',40):
                state['pending'].insert(0,node)
                state['status']='blocked_checkpoint_no_tile_enclosure'
                state['blocked_reason']='maximum subdivision depth reached'
                break
            axis=exc.axis
            if exc.reason=='no finite analytic neighborhood':
                widths=[float(Fraction(hi)-Fraction(lo)) for lo,hi in node['bounds']]
                scale=max(float(abs(origin).upper()),float(abs(origin+delta).upper()),1e-30)
                scales=[float(abs(delta).upper()),scale,scale]
                axis=max(range(3),key=lambda j:widths[j]*scales[j])
            children=split_node(node,axis)
            state['splits'].append({'id':node['id'],'axis':axis})
            state['pending'][0:0]=children
        except (ArithmeticError,ZeroDivisionError,ValueError) as exc:
            state['pending'].insert(0,node)
            state['status']='blocked_checkpoint_no_tile_enclosure'
            state['blocked_reason']=type(exc).__name__+': '+str(exc)
            break
    state['stages'].append({'seconds':time.monotonic()-start,'integrand_calls':calls,'accepted_leaves_added':len(state['accepted'])-accepted_before,'splits_added':len(state['splits'])-splits_before})
    audit_partition(state)
    total=recompute_total(state)
    if not state['pending']:
        state['status']='rigorous_complete_tile_enclosure'
        vals=[x.imag for x in total] if settings['mode']=='qt' else total
        area=(exact(settings['bounds'][1])-exact(settings['bounds'][0]))*(exact(settings['bounds'][3])-exact(settings['bounds'][2]))
        radius=max(float(x.rad()) if settings['mode']=='qt' else max(float(x.real.rad()),float(x.imag.rad())) for x in vals)
        state['moments']=list(map(ball_serialization,vals));state['maximum_radius_diagnostic']=radius
        state['normalized_maximum_radius_diagnostic']=radius/float(area)
        state['claim_limit']='Only the complete specified angular tile; no full-domain or lattice-floor claim.'
    else:
        state['claim_limit']='Accepted leaves only; pending boxes remain. No complete tile enclosure.'
    atomic_json(dest,state)
    return state

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bounds',required=True,help='u0,u1,v0,v1 exact fractions');p.add_argument('--chart',type=int,required=True);p.add_argument('--mode',choices=['qt','green'],required=True);p.add_argument('--height',default='1/4');p.add_argument('--orders',default='5,5,5');p.add_argument('--density-tolerance',default='1/25000');p.add_argument('--bits',type=int,default=128);p.add_argument('--max-seconds',type=float,default=20);p.add_argument('--output',required=True)
    a=p.parse_args()
    settings={'bounds':a.bounds.split(','),'chart':a.chart,'mode':a.mode,'height':a.height,'orders':list(map(int,a.orders.split(','))),'density_tolerance':a.density_tolerance,'bits':a.bits}
    if a.bits<64:p.error('precision must be at least64 bits')
    result=integrate_stage(settings,a.output,a.max_seconds)
    print(json.dumps({'status':result['status'],'accepted':len(result['accepted']),'pending':len(result['pending']),'last_stage':result['stages'][-1],'normalized_radius':result.get('normalized_maximum_radius_diagnostic'),'endpoint':result['endpoint_certificate']['kind']}))
