"""Reconstruct the reference packet and repeat the Section VII analysis of retained data."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

DATA = Path(__file__).resolve().parents[1] / 'evidence' / 'packet_reference'


def read(name):
    return json.loads((DATA / name).read_bytes())


def rows(name):
    with (DATA / name).open() as stream:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]


def check(x, y, tol):
    np.testing.assert_allclose(x, y, rtol=0, atol=tol)


def mean(times, values):
    return float(np.trapezoid(values, times)/(times[-1] - times[0]))


def incident_measure():
    p = read('preparation.json')
    count, a, mass, stiffness = p['n'], p['a'], p['m_q'], p['K_q']
    wave_speed = a * np.sqrt(stiffness / mass)
    ka0 = 2*np.arcsin(p['omega0']*a/(2*wave_speed))
    x = np.arange(count)*a - p['center']
    u = p['amplitude']*np.exp(-x**2/(2*p['sigma']**2))*np.cos(ka0*x/a)
    k = 2*np.pi*np.fft.fftfreq(count, a)
    omega = 2*wave_speed/a * abs(np.sin(k*a/2))
    uh = np.fft.fft(u)
    vh = -1j*np.sign(k)*omega*uh
    vh[[0,count//2]] = 0
    v = np.fft.ifft(vh)
    check(v.imag, 0, 1e-17)
    weights = mass*(abs(vh)**2 + (omega*abs(uh))**2)/(2*count)
    periodic = weights.sum()
    e0 = float(0.5*mass*(v.real @ v.real) + 0.5*stiffness*np.sum(np.diff(u)**2))
    check(periodic-e0, 0.5*stiffness*(u[0]-u[-1])**2, 1e-17)
    # Eight diagonal springs per segment with a bulk initially at rest.
    g = np.array([0.1,0.2,0.2,0.2,0.2,0.1])
    link0 = 2*np.sum(g*(u[60:66]**2 + u[61:67]**2))/3
    check(link0/e0, 2.0502029229148544e-18, 1e-29)
    order = np.argsort(omega)
    cdf = np.cumsum(weights[order])/periodic
    quantiles = omega[order][np.searchsorted(cdf, [0.01,0.99])]
    check(quantiles, [0.9317566688177807,1.4726562160045653], 1e-14)
    return omega, weights/periodic, e0, float(link0/e0), quantiles.tolist(), float(periodic/e0)


def spectral_comparison(omega, weights, direct):
    pos = np.flatnonzero((np.arange(len(omega)) < len(omega)//2)
                         & (omega >= .8) & (omega <= 1.6))
    assert len(pos) == 14 and weights[0] == 0
    w = weights[pos] + weights[-pos]
    coverage = float(w.sum())
    sums = {}
    for key in ('surface_128','surface_160','tetra_96','tetra_128','tetra_160'):
        curve = read(key+'.json')['rows']
        q = []
        for om in omega[pos]:
            pairs = [r for r in curve if abs(r['omega']-om)<1e-12]
            assert len(pairs) == 1
            q.append(pairs[0]['T_qT'])
        sums[key] = float(w @ q)
    check(list(sums.values()), [.1276811,.1276730,.1275676,.1275947,.1276094], 5.1e-8)
    check(coverage, .9996834927581062, 3e-14)
    s = sums['surface_160']
    d = (s-direct)/s*100
    m = (s-sums['tetra_128'])/s*100
    check([d,m], [.0365,.06137], 5.1e-5)
    return {'central_frequencies':len(pos),'periodic_packet_energy_covered':coverage,
            'weighted_in_band_contributions':sums, 'surface_tail_only_range':[s,s+1-coverage],
            'fraction_difference_in_percentage_points':100*(s-direct),
            'prediction_normalized_difference_percent':d,'two_method_difference_percent':m}


def select(geom, modes):
    passing = []
    for d, e in zip(geom, modes, strict=True):
        check(d['time'],e['time'], 1e-11)
        discrepancy = abs(e['periodic_total_over_E0']/d['T_bulk_open'] - 1)
        passing.append(d['T_coupling']<=.0025 and abs(d['current'])<=.0005 and
                       d['historical_shell']<=.002 and discrepancy<=.0025 and abs(d['drift'])<=1e-5)
    starts = [i for i, x in enumerate(passing) if x and (i==0 or not passing[i-1])]
    ends = [i for i, x in enumerate(passing) if x and (i+1==len(passing) or not passing[i+1])]
    candidates = [(geom[i]['time'],geom[j]['time']) for i,j in zip(starts,ends)
                  if geom[j]['time']-geom[i]['time']>=2-1e-12]
    return max(candidates, key=lambda ij:(round(ij[1]-ij[0],8),-ij[0]))


def summary(modal, interval, spacing):
    t0, t1 = interval
    selected = [r for r in modal if t0-1e-8 <= r['time'] <= t1+1e-8
                and abs((r['time']-t0)/spacing-round((r['time']-t0)/spacing))<1e-8]
    t = [r['time'] for r in selected]
    check([t[0], t[-1]], [t0,t1], 1e-10)
    low = mean(t,[r['qT_min_over_E0'] for r in selected])
    high = mean(t,[r['qT_max_over_E0'] for r in selected])
    total = mean(t,[r['periodic_total_over_E0'] for r in selected])
    for r in selected:
        check(r['periodic_total_over_E0'],r['qT_min_over_E0']+r['qL_max_over_E0'],5e-14)
    return low,high,(low+high)/2,total


def dynamics(link0_fraction):
    m = read('REFERENCE_MODAL_ANALYSIS.json')['rows']
    d = rows('REFERENCE_TRAJECTORY_ROWS.csv')
    assert len(m)==len(d)==46
    window = select(d,m)
    check(window, [52.8,55.4],1e-12)
    low,high,direct,total=summary(m,window,.2)
    check([low,high,total],[.12762638643388616,.12762638738500592,.13892844304455843],5e-14)
    dsel = [r for r in d if window[0]-1e-9<=r['time']<=window[1]+1e-9]
    check(len(dsel),14,0)
    t=np.array([r['time'] for r in dsel])
    inc=np.array([r['T_current'] for r in dsel])
    gain=np.array([r['T_bulk_open']+(r['T_coupling']-link0_fraction)/2 for r in dsel])
    max_residual=float(100*np.max(abs(inc/gain-1)))
    check(max_residual,.0296799,5e-8)
    # Independently integrate every stored integrator-level instantaneous current in
    # the reporting interval; the first cumulative value includes the earlier drive.
    full=rows('FULL_STEP_DIRECT_THRESHOLDS.csv')
    tt=np.array([r['time'] for r in full]); jj=np.array([r['current'] for r in full])
    cumulative=inc[0]+np.r_[0,np.cumsum(np.diff(tt)*(jj[1:]+jj[:-1])/2)]
    check(cumulative,[r['T_current'] for r in full],2e-15)
    check(len(full),66,0)
    drift=max(abs(r['drift']) for r in full)
    check(drift,3.54e-7,5.1e-10)
    shell=max(r['union'] for r in read('REFERENCE_SHELL_AUDIT.json')['rows'])
    check(shell,.001858,5.1e-7)
    small_step=select(d[::2],m[::2])
    check(small_step,[52.8,55.2],2e-12)
    coarse=summary(m,small_step,.4)[2]
    fine_common=summary(m,small_step,.2)[2]
    common={variant:summary(data,(52.8,55.2),.6)[2] for variant,data in
            [('0.03',read('DT003_MODAL_CORRECTED.json')['rows']),('0.04',m),('0.05',read('DT005_MODAL_CORRECTED.json')['rows'])]}
    span=(max(common.values())-min(common.values()))/common['0.04']*100
    check(span,.00160,5e-6)
    wide=summary(read('BOX64_MODAL_ANALYSIS.json')['rows'],window,.2)
    box=100*abs(wide[2]/direct-1)
    ql_box=100*abs((wide[3]-wide[2])/(total-direct)-1)
    check([box,ql_box],[.04398,.7938],.00005)
    wide_rows=[r for r in rows('BOX64_ROWS.csv') if window[0]-1e-9<=r['time']<=window[1]+1e-9]
    avg_box={key:abs(mean(t,[r[key] for r in wide_rows])/mean(t,[r[key] for r in dsel])-1)
             for key in ['T_bulk_open','T_current']}
    assert max(avg_box.values()) < 4.4e-8
    return {'reporting_window':list(window),'samples_in_window':len(dsel),
            'qT_allocation_interval':[low,high],'direct_qT_fraction':direct,
            'qT_fraction_of_periodic_bulk':direct/total,'periodic_bulk_fraction':total,
            'max_current_relative_residual_percent':max_residual,'max_fullstep_energy_drift_fraction':drift,
            'postselection_total_six_face_shell_fraction':shell,
            'coarsened_window':list(small_step),'coarse_window_change_percent':100*abs(coarse/direct-1),
            'common_interval_cadence_change_percent':100*abs(coarse/fine_common-1),
            'three_step_sizes_qT_relative_span_percent':span,
            'fixed_window_larger_box_qT_change_percent':box,'fixed_window_larger_box_qL_change_percent':ql_box,
            'larger_box_open_energy_and_current_relative_changes':avg_box}


def main():
    for name, item in read('PROVENANCE.json')['files'].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==item['sha256'],name
    om,w,e0,link0,quantiles,ratio=incident_measure()
    trajectory=dynamics(link0)
    spectral=spectral_comparison(om,w,trajectory['direct_qT_fraction'])
    fig=read('figure6a_data.json')
    m=read('REFERENCE_MODAL_ANALYSIS.json')['rows']
    check(fig['time'], [r['time'] for r in m],1e-12)
    check(fig['qT'], [r['qT_production_over_E0'] for r in m],2e-15)
    check(fig['accepted_interval'],trajectory['reporting_window'],1e-12)
    check(fig['prediction_range'],[spectral['weighted_in_band_contributions'][k] for k in ['tetra_128','surface_160']],2e-15)
    print(json.dumps({'status':'PASS','preparation':{'initial_open_line_energy':e0,
        'initial_link_energy_over_open_line_energy':link0,'periodic_to_open_energy_ratio':ratio,
        'Nyquist_weight':float(w[len(w)//2]),'one_to_99_percent_omega':quantiles},'trajectory':trajectory,
        'stationary':spectral,'figure_6a_original_curve_points_matched':len(m),
        'scope':'DFT preparation reconstructed; reference calculations reanalysed from original step and snapshot reductions. Optional trajectory reproduction and independent Fourier projection have a separate receipt.'},indent=2))


if __name__=='__main__':
    main()
