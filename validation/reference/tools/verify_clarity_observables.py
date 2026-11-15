"""Check the energy denominators, qT octant result and contact comparison in Section VI."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def read(directory, filename):
    return json.loads((directory / filename).read_text())


def octant_result(raw, reference, grid):
    """Compute the qT distribution directly from the two lower rows of the 3-by-8 matrix."""
    assert raw['grid'] == grid and raw['omega'] == 1.2
    powers = raw['branch_octant_powers']
    assert len(powers) == 3 and all(len(row) == 8 for row in powers)
    assert all(power > 0 for row in powers for power in row)
    octant_qt = [p0 + p1 for p0, p1 in zip(*powers[:2])]
    qt = sum(octant_qt)
    total = sum(map(sum, powers))
    assert math.isclose(qt, raw['reconstructed']['T_qT'], rel_tol=1e-12)
    assert math.isclose(total, raw['reconstructed']['Tbulk'], rel_tol=1e-12)
    smallest = 100 * min(octant_qt) / qt
    smallest_bulk = 100 * min(map(sum, zip(*powers))) / total
    assert round(smallest, 1) == 1.4 and round(smallest_bulk, 1) == 1.3
    assert abs(total / reference['Tbulk'] - 1) < 0.00022
    assert round(100 * raw['off_axis']['fraction_sin_theta_gt_0p5']) == 94
    assert round(raw['off_axis']['weighted_median_sin_theta'], 2) == 0.86
    return {'least_octant_percent_of_emitted_qT':smallest, 'least_octant_percent_of_total_bulk':smallest_bulk}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=Path(__file__).resolve().parents[1]/'evidence')
    evidence=parser.parse_args().evidence
    data=evidence/'clarity'
    for f, record in read(data,'PROVENANCE.json')['files'].items():
        assert hashlib.sha256((data/f).read_bytes()).hexdigest() == record['sha256']
    source=read(data,'DOWNSTREAM_SURFACE_SELECTED.json')
    clusters=source['reference']
    geo=source['same_state']['wave_normal_geometric']['total']
    assert round(geo,7) == 0.1256316
    assert round(100*geo,1) == round(100*clusters['TT_modal'],1) == 12.6
    assert abs(geo-clusters['TT_modal']) > 0.00002
    octants={grid:octant_result(read(data,f'ANGULAR_SURFACE_SELECTED_N{grid}.json'),clusters,grid) for grid in (80,96)}

    fs=evidence/'vi_results'
    s=read(fs,'SCALE_PROBE.json')['states']
    ramp=read(fs,'PROFILE_COMPARISON.json')['best_sampled_refined']
    controls=[next(row for row in s if row['cache']=='N192_M256' and row['scale']==v) for v in (1.,2.)]
    rows=controls+[ramp]
    power_partition=[]
    for r in rows:
        qt, ql=r['TqT'],r['TqL']
        bulk=qt+ql
        assert math.isclose(bulk,r['Tbulk'],abs_tol=1e-12)
        assert abs(1-bulk-r['Tline']-r['Rq'])<2e-11
        power_partition.append({'percent_incident_emitted':100*bulk,'qT_percent_emitted':100*qt/bulk,'qT_percent_incident':100*qt})
    assert [round(x['qT_percent_incident'],1) for x in power_partition]==[12.6,52.2,73.6]
    assert [round(x['percent_incident_emitted'],1) for x in power_partition]==[13.7,57.2,80.4]
    assert [round(x['qT_percent_emitted'],1) for x in power_partition]==[91.9,91.3,91.5]
    cones=read(evidence/'controls','integral_fine.json')
    full=cones['cone_branch_powers']['group_forward_10']['full']
    all_qt=sum(cones['total_branch_powers']['full'][:2])
    qt_cone=sum(full[:2])
    frac_qt_of_cone=qt_cone/sum(full)
    capture=qt_cone/all_qt
    assert round(100*frac_qt_of_cone,2)==1.19
    assert round(100*capture,4)==0.0073
    print(json.dumps({'status':'PASS','octants':octants,
                      'contact_order':['reference taper','doubled taper','increasing ramp'],
                      'power_partition':power_partition,
                      'ten_degree_forward_cone':{'incident_fraction_qT':qt_cone,'qT_percent_collected':100*frac_qt_of_cone,'percent_of_all_emitted_qT_captured':100*capture},
                      'same_state_wave_normal_projection':geo,'same_state_modal_qT':clusters['TT_modal'],
                      'scope':'Recomputed ratios and branch/octant sums from stored integrations; angular data are numerical quadratures.'},indent=2))


if __name__=='__main__':
    main()
