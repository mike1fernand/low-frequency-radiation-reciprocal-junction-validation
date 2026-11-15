"""Regenerate the open 120 x 80 x 80 ramp with its archived C++17/NumPy producers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from verify_ramp_packet import EVIDENCE, read, require, verify_archive_inputs


def digest(path):
    checksum = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024*1024), b''):
            checksum.update(chunk)
    return checksum.hexdigest()


def run_producer(output):
    # The reproduced configuration is the 80-box primary. The historical 56-box
    # trial and 96-box/time-step controls are checked by verify_ramp_packet.py.
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
    engine = output / 'ramp_engine'
    data_dir = output / 'packet' / 'candidate_dt004_box80'
    data_dir.mkdir(parents=True)
    review = output / 'packet_review'
    review.mkdir()
    subprocess.run(['g++', '-O3', '-std=c++17', '-fopenmp', EVIDENCE/'evolve_candidate.cpp',
                    '-o', engine], check=True)
    cfg = read('configuration.json')
    command = [str(engine), str(data_dir), str(cfg['target_final_step']), str(cfg['dt']),
               str(cfg['quarter_shape'][1]), *[format(g,'.17g') for g in cfg['profile']], '1650']
    with (EVIDENCE/'initial_qv.txt').open('r') as preparation:
        run = subprocess.run(command, stdin=preparation, capture_output=True, text=True, env=env, check=True)
    evolution = json.loads(run.stdout)
    require(evolution['status'] == 'PASS' and evolution['start_step'] == 0 and evolution['target_step'] == 1650,
            'Incomplete ramp evolution')
    shutil.copyfile(EVIDENCE/'configuration.json', data_dir/'configuration.json')
    analysis = subprocess.check_output([sys.executable, EVIDENCE/'analyze_candidate.py', data_dir], env=env)
    (output/'fresh_modal_analysis.json').write_bytes(analysis)
    actual = json.loads(analysis)
    archived = read('CANDIDATE_PACKET_BOX80_RESULT.json')
    require(actual['longest'] == archived['longest'] and actual['passing_runs'] == archived['passing_runs'],
            'The five-condition selector differs from the archived result')
    for a,b in zip(actual['rows'],archived['rows'],strict=True):
        require(a['passed']==b['passed'] and a['t']==b['t'], 'Selector rows differ')
        for key in ('open_bulk','link','J','shell','periodic_difference','drift'):
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-11)
    np.testing.assert_allclose(actual['branch_samples'], archived['branch_samples'], atol=1e-11, rtol=0)
    shutil.copyfile(EVIDENCE/'CANDIDATE_PACKET_BOX80_RESULT.json', output/'packet'/'CANDIDATE_PACKET_BOX80_RESULT.json')
    shutil.copyfile(EVIDENCE/'check_box80_analysis.py', review/'check_box80_analysis.py')
    independent = json.loads(subprocess.check_output([sys.executable,review/'check_box80_analysis.py'], env=env))
    return actual, independent, data_dir


def receipt(result, independent, directory):
    raw = read('RAW_DELIVERY_SCOPE.json')
    hashes = {}
    for name in ['line.bin','corner.bin','snap_times.bin','snap.bin','checkpoint_0.bin','checkpoint_1650.bin']:
        record = next(row for row in raw['files'] if row['path']=='candidate_dt004_box80/'+name)
        hashes[name] = (digest(directory/name) == record['sha256'])
    candidates = independent['three_snapshot_checks']
    return {'status':'PASS','evolution':'Fresh full run with the archived producer and input.',
            'physical_bulk_shape':[120,80,80], 'time_step':0.04, 'steps':1650,
            'archived_raw_hash_matches': hashes,
            'selector_rows_recalculated':len(result['rows']), 'window':result['longest'],
            'qT_allocations_recalculated':len(result['branch_samples']),
            'mean_qT_lower_upper':result['mean_low_high_total'][:2],
            'qT_midpoint':sum(result['mean_low_high_total'][:2])/2,
            'independent_snapshot_times':[r['time'] for r in candidates],
            'max_independent_snapshot_difference':max(r['maximum_abs_difference'] for r in candidates),
            'scope':'The original dynamical and analysis codes were re-executed. A separate spectral implementation checked three snapshots. Other box/time-step controls are archival in this revision.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True,help='New directory for binary histories (about 1 GiB).')
    args=parser.parse_args()
    verify_archive_inputs()
    result, independent, directory = run_producer(args.output_dir.resolve())
    got = receipt(result, independent, directory)
    print(json.dumps(got,indent=2))


if __name__=='__main__':
    main()
