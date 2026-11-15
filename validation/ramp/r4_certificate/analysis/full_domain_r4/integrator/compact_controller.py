"""Bounded restartable rigorous integration, with exact compact stage sums.

The three copied R3 primitive files are unchanged. Vector tolerances use
positive constant scaling before that primitive and rescaling afterwards.
Exact endpoint addition prevents widening caused by repeated serialization.
"""
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import arb, acb, ctx
from arbitrary_tile_sharp import endpoint_certificate, contour_legs, split_node
from certified_tensor_sharp import gauss_box, NeighborhoodFailure
from source_moments import moment_values, ball_serialization

SCHEMA = 'paper1-compact-tile-v1'
SOURCE_NAMES = ('compact_controller.py', 'arbitrary_tile_sharp.py',
                'certified_tensor_sharp.py', 'source_moments.py')


def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(x):
    return hashlib.sha256(canonical(x)).hexdigest()


def exact(x):
    q = Fraction(x)
    return arb(q.numerator) / q.denominator


def normalize_endpoint(m, e):
    m, e = int(m), int(e)
    if not m:
        return ['0', 0]
    shifts = (abs(m) & -abs(m)).bit_length() - 1
    return [str(m >> shifts), e + shifts]


def add_endpoint(a, b):
    am, ae = int(a[0]), int(a[1])
    bm, be = int(b[0]), int(b[1])
    e = min(ae, be)
    return normalize_endpoint((am << (ae-e)) + (bm << (be-e)), e)


def endpoint_arb(a):
    return arb(int(a[0])) * arb(2) ** int(a[1])


def endpoint_fraction(a):
    m, e = int(a[0]), int(a[1])
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)


def zero_sum():
    # Each moment: real lower, real upper, imaginary lower, imaginary upper.
    return [[['0', 0] for _ in range(4)] for _ in range(17)]


def ball_endpoints(vals):
    if len(vals) != 17 or not all(x.is_finite() for x in vals):
        raise ValueError('exactly 17 finite complex values required')
    return [[normalize_endpoint(*y.man_exp()) for y in
             (x.real.lower(), x.real.upper(), x.imag.lower(), x.imag.upper())]
            for x in vals]


def add_sums(a, b):
    if len(a) != 17 or len(b) != 17:
        raise ValueError('wrong moment count')
    return [[add_endpoint(x, y) for x, y in zip(ar, br)]
            for ar, br in zip(a, b)]


def validate_sum(values):
    if len(values) != 17:
        raise ValueError('wrong aggregate moment count')
    for row in values:
        if len(row) != 4:
            raise ValueError('wrong aggregate endpoint count')
        for end in row:
            if len(end) != 2 or end != normalize_endpoint(*end):
                raise ValueError('noncanonical dyadic endpoint')
        if (endpoint_fraction(row[0]) > endpoint_fraction(row[1]) or
                endpoint_fraction(row[2]) > endpoint_fraction(row[3])):
            raise ValueError('reversed interval endpoints')


def aggregate_balls(values):
    validate_sum(values)
    return [acb(endpoint_arb(x[0]).union(endpoint_arb(x[1])),
                endpoint_arb(x[2]).union(endpoint_arb(x[3]))) for x in values]


def normalized_settings(settings):
    s = dict(settings)
    allowed = {'bounds', 'chart', 'mode', 'height', 'orders', 'bits',
               'density_tolerance', 'density_tolerances', 'max_depth',
               'radius_multipliers', 'endpoint_candidates'}
    if set(s) - allowed:
        raise ValueError('unknown settings: ' + repr(sorted(set(s)-allowed)))
    s['bounds'] = [str(Fraction(x)) for x in s['bounds']]
    if len(s['bounds']) != 4:
        raise ValueError('four angular bounds required')
    if not isinstance(s.get('bits', 128), int) or s.get('bits', 128) < 64:
        raise ValueError('precision must be an integer >=64')
    s['bits'] = s.get('bits', 128)
    s['orders'] = list(s.get('orders', [5, 5, 5]))
    if len(s['orders']) != 3 or any(x not in (3, 5) for x in s['orders']):
        raise ValueError('three orders, each 3 or 5, required')
    if 'density_tolerance' in s and 'density_tolerances' in s:
        raise ValueError('choose scalar or vector tolerance')
    tolerances = s.pop('density_tolerance', None)
    if tolerances is not None:
        tolerances = [tolerances] * 17
    else:
        tolerances = s.get('density_tolerances', ['1/25000'] * 17)
    if len(tolerances) != 17:
        raise ValueError('17 density tolerances required')
    s['density_tolerances'] = [str(Fraction(x)) for x in tolerances]
    if any(Fraction(x) <= 0 for x in s['density_tolerances']):
        raise ValueError('positive density tolerances required')
    s['height'] = str(Fraction(s.get('height', '1/4')))
    s['max_depth'] = s.get('max_depth', 40)
    if not isinstance(s['max_depth'], int) or s['max_depth'] < 1:
        raise ValueError('positive integer max_depth required')
    s['radius_multipliers'] = list(s.get('radius_multipliers', [2, 4, 8, 16]))
    if not s['radius_multipliers'] or any(not isinstance(x, int) or x <= 1 for x in s['radius_multipliers']):
        raise ValueError('radius multipliers must be integers >1')
    s['endpoint_candidates'] = [str(Fraction(x)) for x in s.get('endpoint_candidates', ['1/2','9/16','5/8','11/16'])]
    return s


def identity():
    return {'sources': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                        for name in SOURCE_NAMES},
            'python': sys.version, 'flint_version': __import__('flint').__version__}


def geometry(state):
    s = state['settings']
    expected_legs = {'qt': 3, 'green': 4}.get(s['mode'])
    if state['leg_count'] != expected_legs:
        raise ValueError('contour leg count inconsistent with mode')
    a, b, c, d = s['bounds']
    nodes = {f'{j}:': {'id':f'{j}:', 'leg':j,
                      'bounds':[['0','1'],[a,b],[c,d]], 'depth':0}
             for j in range(expected_legs)}
    for node_id, axis in state['splits']:
        if axis not in (0,1,2) or node_id not in nodes:
            raise ValueError('invalid split tree')
        node = nodes.pop(node_id)
        nodes.update({x['id']:x for x in split_node(node,axis)})
    leaves = state['accepted_ids'] + state['pending_ids']
    if len(leaves) != len(set(leaves)) or set(leaves) != set(nodes):
        raise ValueError('missing, duplicate or unexpected coverage leaf')
    return nodes


def audit_state(state):
    if state['schema'] != SCHEMA:
        raise ValueError('unknown checkpoint schema')
    if state.get('state_sha256') != digest({k:v for k,v in state.items() if k != 'state_sha256'}):
        raise ValueError('checkpoint integrity mismatch')
    nodes = geometry(state)
    total, previous, position = zero_sum(), None, 0
    for index, stage in enumerate(state['stages']):
        payload = {k:v for k,v in stage.items() if k != 'sha256'}
        if stage['sha256'] != digest(payload) or stage['previous_stage_sha256'] != previous:
            raise ValueError('stage hash chain mismatch')
        if stage['index'] != index or stage['aggregate_before_sha256'] != digest(total):
            raise ValueError('stage aggregate chain mismatch')
        count = stage['accepted_count']
        if not isinstance(count,int) or count < 0:
            raise ValueError('invalid accepted stage count')
        accepted_ids = state['accepted_ids'][position:position+count]
        if len(accepted_ids) != count or digest(accepted_ids) != stage['accepted_ids_sha256']:
            raise ValueError('stage accepted membership mismatch')
        validate_sum(stage['increment_endpoints'])
        total = add_sums(total, stage['increment_endpoints'])
        if stage['aggregate_after_sha256'] != digest(total):
            raise ValueError('stage aggregate mismatch')
        position += count
        previous = stage['sha256']
    validate_sum(state['aggregate_endpoints'])
    if position != len(state['accepted_ids']) or total != state['aggregate_endpoints']:
        raise ValueError('compact aggregate disagrees with stage increments')
    complete = not state['pending_ids']
    if (state['status'] == 'rigorous_complete_tile_enclosure') != complete:
        raise ValueError('completion status inconsistent with coverage')
    if complete:
        # Audit serialized public moments independently of their display fields.
        ctx.prec = state['settings']['bits']
        vals = aggregate_balls(total)
        if state['settings']['mode'] == 'qt':
            vals = [x.imag for x in vals]
        if state.get('moments') != list(map(ball_serialization, vals)):
            raise ValueError('public moments disagree with compact exact aggregate')
        for key, value in exact_radius_diagnostics(state).items():
            if state.get(key) != value:
                raise ValueError('radius diagnostics disagree with exact aggregate')
    return nodes


def atomic_json(destination, state):
    state.pop('state_sha256',None)
    state['state_sha256'] = digest(state)
    destination = Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    tmp = destination.with_suffix(destination.suffix+'.partial')
    with tmp.open('w') as stream:
        json.dump(state,stream,separators=(',',':'),allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp,destination)


def exact_radius_diagnostics(state):
    area = ((Fraction(state['settings']['bounds'][1])-Fraction(state['settings']['bounds'][0])) *
            (Fraction(state['settings']['bounds'][3])-Fraction(state['settings']['bounds'][2])))
    mode = state['settings']['mode']
    radii = []
    for row in state['aggregate_endpoints']:
        rre = (endpoint_fraction(row[1])-endpoint_fraction(row[0])) / 2
        rim = (endpoint_fraction(row[3])-endpoint_fraction(row[2])) / 2
        radii.append(rim if mode == 'qt' else max(rre,rim))
    return {'radii_diagnostic':list(map(float,radii)),
            'normalized_radii_diagnostic':[float(x/area) for x in radii],
            'maximum_radius_diagnostic':float(max(radii)),
            'normalized_maximum_radius_diagnostic':float(max(radii)/area)}


def completed_moments(state):
    audit_state(state)
    if state['pending_ids']:
        raise ValueError('complete claim has pending boxes')
    total = aggregate_balls(state['aggregate_endpoints'])
    return [x.imag for x in total] if state['settings']['mode']=='qt' else total


def integrate_stage(settings, destination, max_seconds=20):
    if not math.isfinite(max_seconds) or not 0 < max_seconds <= 40:
        raise ValueError('stage budget must be positive and <=40 seconds')
    settings = normalized_settings(settings)
    ctx.prec = settings['bits']
    t0, cert = endpoint_certificate(settings['bounds'],settings['chart'],settings['endpoint_candidates'])
    legs = contour_legs(settings['mode'],settings['height'],t0)
    dest = Path(destination)
    current_identity = identity()
    if dest.exists():
        state = json.loads(dest.read_text())
        if state['settings'] != settings or state['identity'] != current_identity:
            raise ValueError('checkpoint source/runtime/settings mismatch')
        if state['endpoint_certificate'] != cert:
            raise ValueError('endpoint certificate differs from recomputation')
        nodes = audit_state(state)
        if state['status']=='rigorous_complete_tile_enclosure':
            return state
    else:
        state = {'schema':SCHEMA,'status':'partial_checkpoint_no_tile_enclosure',
                 'settings':settings,'identity':current_identity,
                 'endpoint_certificate':cert,'leg_count':len(legs),
                 'pending_ids':[f'{j}:' for j in range(len(legs))],
                 'accepted_ids':[],'splits':[],'aggregate_endpoints':zero_sum(),'stages':[]}
        nodes = geometry(state)
    start = time.monotonic()
    calls = 0
    accepted_before = len(state['accepted_ids'])
    splits_before = len(state['splits'])
    before = state['aggregate_endpoints']
    increment = zero_sum()
    scales = list(map(exact,settings['density_tolerances']))
    state['status'] = 'partial_checkpoint_no_tile_enclosure'
    state.pop('blocked_reason',None)
    while state['pending_ids'] and time.monotonic()-start < max_seconds:
        node_id = state['pending_ids'].pop(0)
        node = nodes[node_id]
        _,origin,delta = legs[node['leg']]
        def callback(q,u,v):
            nonlocal calls
            if time.monotonic()-start >= max_seconds:
                raise TimeoutError('stage time limit')
            calls += 1
            return [x*delta/scales[m] for m,x in enumerate(moment_values(origin+delta*q,u,v,settings['chart']))]
        try:
            vals, meta = gauss_box(callback,[(exact(lo),exact(hi)) for lo,hi in node['bounds']],
                                  tuple(settings['orders']),tuple(settings['radius_multipliers']),arb(1)/len(legs))
            vals = [x*scales[m] for m,x in enumerate(vals)]
            increment = add_sums(increment,ball_endpoints(vals))
            state['accepted_ids'].append(node_id)
        except TimeoutError:
            state['pending_ids'].insert(0,node_id)
            break
        except NeighborhoodFailure as exc:
            if node['depth'] >= settings['max_depth']:
                state['pending_ids'].insert(0,node_id)
                state['status'] = 'blocked_checkpoint_no_tile_enclosure'
                state['blocked_reason'] = 'maximum subdivision depth reached'
                break
            axis = exc.axis
            if exc.reason == 'no finite analytic neighborhood':
                widths = [float(Fraction(hi)-Fraction(lo)) for lo,hi in node['bounds']]
                scale = max(float(abs(origin).upper()),float(abs(origin+delta).upper()),1e-30)
                axis = max(range(3),key=lambda j:widths[j]*[float(abs(delta).upper()),scale,scale][j])
            children = split_node(node,axis)
            del nodes[node_id]
            nodes.update({x['id']:x for x in children})
            state['splits'].append([node_id,axis])
            state['pending_ids'][0:0] = [x['id'] for x in children]
        except (ArithmeticError,ZeroDivisionError,ValueError) as exc:
            state['pending_ids'].insert(0,node_id)
            state['status'] = 'blocked_checkpoint_no_tile_enclosure'
            state['blocked_reason'] = type(exc).__name__+': '+str(exc)
            break
    state['aggregate_endpoints'] = add_sums(before,increment)
    added_ids = state['accepted_ids'][accepted_before:]
    stage = {'index':len(state['stages']),
             'previous_stage_sha256':state['stages'][-1]['sha256'] if state['stages'] else None,
             'seconds':time.monotonic()-start,'integrand_calls':calls,
             'accepted_count':len(added_ids),'accepted_ids_sha256':digest(added_ids),
             'splits_added':len(state['splits'])-splits_before,
             'aggregate_before_sha256':digest(before),'increment_endpoints':increment,
             'aggregate_after_sha256':digest(state['aggregate_endpoints'])}
    stage['sha256'] = digest(stage)
    state['stages'].append(stage)
    if not state['pending_ids']:
        state['status'] = 'rigorous_complete_tile_enclosure'
        vals = aggregate_balls(state['aggregate_endpoints'])
        if settings['mode']=='qt':
            vals = [x.imag for x in vals]
        state['moments'] = list(map(ball_serialization,vals))
        state.update(exact_radius_diagnostics(state))
        state['claim_limit'] = 'Complete specified angular tile only. No global or lattice-floor claim.'
    else:
        state['claim_limit'] = 'Accepted leaves only; pending leaves remain. No complete tile enclosure.'
    state.pop('state_sha256',None)
    state['state_sha256'] = digest(state)
    audit_state(state)
    atomic_json(dest,state)
    return state


def summary(state):
    return {'status':state['status'],'accepted':len(state['accepted_ids']),
            'pending':len(state['pending_ids']),'stages':len(state['stages']),
            'last_stage':{k:v for k,v in state['stages'][-1].items() if k not in ('increment_endpoints',)},
            'normalized_radii_diagnostic':state.get('normalized_radii_diagnostic'),
            'endpoint':state['endpoint_certificate']['kind']}
