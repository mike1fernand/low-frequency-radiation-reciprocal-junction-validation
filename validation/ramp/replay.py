"""Run declared validation cases in NEW external worktrees; never edit evidence.

Use one of the separate environments documented in README.md. This runner does
not install dependencies or change scientific parameters. Saved-byte checks do
not replace proof or numerical re-execution.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('task', choices=['response', 'baseline', 'spectrum',
                                   'quadrature', 'packet', 'resume-tests', 'report'])
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--mode', choices=['pilot', 'refine', 'independent-tetra',
                                     'fresh-tetra-n48'], default='pilot')
    ap.add_argument('--case', choices=['base', 'all'], default='base')
    args = ap.parse_args()
    out = args.output_dir.resolve()
    require(not out.exists(), 'Choose a NEW output directory; existing outputs are never overwritten.')
    require(out != ROOT and ROOT not in out.parents,
            'Output must be outside the frozen validation repository.')
    env = dict(os.environ, PYTHONOPTIMIZE='0', PYTHONDONTWRITEBYTECODE='1',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2', MKL_NUM_THREADS='1',
               MPLCONFIGDIR=str(out / 'matplotlib-cache'))
    # Python -B and a cleared optimization flag protect assertion guards in
    # historical source, even when this outer runner was invoked with -O.
    python = [sys.executable, '-B']
    verifier = python + [str(ROOT / 'verify_repository.py')]
    pre = subprocess.run(verifier, env=env, capture_output=True, text=True)
    require(pre.returncode == 0, 'Pre-replay repository check failed:\n' + pre.stdout + pre.stderr)
    out.mkdir(parents=True)
    (out / 'PRECHECK.log').write_text(pre.stdout + pre.stderr)
    commands = []
    start = time.perf_counter()

    def run(argv, label, cwd=None):
        tick = time.perf_counter()
        result = subprocess.run([str(x) for x in argv], cwd=cwd, env=env,
                                capture_output=True, text=True)
        (out / (label + '.stdout')).write_text(result.stdout)
        (out / (label + '.stderr')).write_text(result.stderr)
        commands.append(dict(label=label, argv=[str(x) for x in argv],
                             returncode=result.returncode,
                             wall_seconds=time.perf_counter() - tick))
        require(result.returncode == 0, 'Replay failed: ' + label + '; see external logs in ' + str(out))
        return result

    task_result = {}
    error = None
    completed = False
    try:
        if args.task == 'response':
            run(python + [ROOT / 'tolerance/reproduce_response_certificates.py',
                          '--output-dir', out / 'response'], 'response')
        elif args.task == 'baseline':
            work = out / 'r4_certificate'
            shutil.copytree(ROOT / 'r4_certificate', work)
            ledger = work / 'analysis/full_domain_r4/budget/FINAL_CAMPAIGN_LEDGER.json'
            run(python + [work / 'analysis/full_domain_r4/export/verify_saved.py',
                          '--output', out / 'baseline-result',
                          '--expected-ledger', ledger, '--physical-check', '--bridge'],
                'baseline', cwd=work)
            task_result['scope'] = ('Saved certificate states, fresh compact physical checks and bridge; '
                                    'not reintegration of all 449,698 leaves.')
        elif args.task == 'spectrum':
            run(python + [ROOT / 'packet_spectrum/reproduce.py', '--output-dir',
                          out / 'spectrum'], 'spectrum')
        elif args.task == 'quadrature':
            run(python + [ROOT / 'fresh_kernel/reproduce.py', '--mode', args.mode,
                          '--output', out / 'quadrature'], 'quadrature')
        elif args.task == 'resume-tests':
            run(python + [ROOT / 'packet_guard/test_safe_resume.py', '--packet-dir',
                          ROOT / 'packet', '--output', out / 'guard-tests'], 'resume-tests')
        elif args.task == 'packet':
            work = out / 'packet'
            work.mkdir()
            for name in ['evolve_candidate.cpp', 'run_control.py', 'analyze_candidate.py',
                         'initial_qv.txt', 'PREDECLARED_PACKET_CONFIGURATION.json']:
                shutil.copyfile(ROOT / 'packet' / name, work / name)
            run(['g++', '-std=c++17', '-O3', '-fopenmp', work / 'evolve_candidate.cpp',
                 '-o', work / 'evolve_candidate'], 'packet-compile')
            cases = [('candidate_dt004_box80', 40, '.04', '66',
                      'CANDIDATE_PACKET_BOX80_RESULT.json', False)]
            if args.case == 'all':
                cases += [('control_dt002_box80', 40, '.02', '64.2',
                           'CONTROL_DT002_BOX80_RESULT.json', True),
                          ('control_dt004_box96', 48, '.04', '64.2',
                           'CONTROL_DT004_BOX96_RESULT.json', True)]
            comparisons = []
            for name, half_box, dt, end, frozen_result, control in cases:
                argv = python + [work / 'run_control.py', '--half-box', str(half_box),
                                 '--dt', dt, '--end', end, '--name', name]
                if control:
                    argv += ['--frozen-window', '57.4', '64.2']
                run(argv, name + '-evolution', cwd=work)
                analysis = run(python + [work / 'analyze_candidate.py', work / name],
                               name + '-analysis', cwd=work)
                new = json.loads(analysis.stdout)
                old = json.loads((ROOT / 'packet' / frozen_result).read_text())
                require(new['status'] == 'PASS', 'Packet analysis did not pass: ' + name)
                require(new.get('frozen_window_all_predicates_pass') is True,
                        'A selected-window packet predicate failed: ' + name)
                require(new['longest'] == old['longest'] == [57.4, 64.2],
                        'Packet accepted interval changed: ' + name)
                delta = max(abs(float(a) - float(b)) for a, b in
                            zip(new['mean_low_high_total'], old['mean_low_high_total']))
                require(delta <= 1e-9, 'Packet arithmetic regression differs by >1e-9: ' + name)
                (work / frozen_result).write_text(analysis.stdout)
                comparisons.append(dict(case=name, status='PASS', max_mean_difference=delta,
                                        arithmetic_regression_tolerance=1e-9,
                                        tolerance_scope='Replay consistency; not a physical error enclosure.'))
            task_result['comparisons'] = comparisons
            task_result['cadence_note'] = ('The current analyzer includes the final cadence endpoint; '
                                         'legacy secondary common-mean fields are superseded.')
        elif args.task == 'report':
            work = out / 'report'
            work.mkdir()
            names = ['tolerance/VALIDATED_PROOF_CHAIN.json', 'fresh_kernel/SUMMARY.json',
                     'packet_spectrum/PACKET_SPECTRUM.json',
                     'packet/CANDIDATE_PACKET_BOX80_RESULT.json',
                     'packet/CONTROL_DT002_BOX80_RESULT.json',
                     'packet/CONTROL_DT004_BOX96_RESULT.json',
                     'assets/fonts/DejaVuSans.ttf', 'assets/fonts/DejaVuSans-Bold.ttf']
            for name in names:
                target = work / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / name, target)
            source = ROOT / 'make_validation_report.py'
            text = source.read_text()
            substitutions = []
            for font in ['DejaVuSans.ttf', 'DejaVuSans-Bold.ttf']:
                old = "'/usr/share/fonts/truetype/dejavu/" + font + "'"
                new = "str(T/'assets/fonts/" + font + "')"
                require(text.count(old) == 1, 'Report font path anchor changed: ' + font)
                text = text.replace(old, new, 1)
                substitutions.append(dict(old=old, new=new))
            staged = work / source.name
            staged.write_text(text)
            task_result['io_adapter'] = dict(original_sha256=digest(source),
                                             staged_sha256=digest(staged),
                                             exact_substitutions=substitutions)
            run(python + [staged], 'report', cwd=work)
            task_result['scope'] = ('Rebuild of the dated four-page scientific report; its historical '
                                    'split-delivery paragraph is superseded by the current root README.')
        completed = True
    except BaseException as exc:
        error = type(exc).__name__ + ': ' + str(exc)
    finally:
        if not completed and error is None:
            error = 'Task did not complete.'
        post = subprocess.run(verifier, env=env, capture_output=True, text=True)
        (out / 'POSTCHECK.log').write_text(post.stdout + post.stderr)
        if post.returncode:
            error = (error + '\n' if error else '') + 'Post-replay frozen-repository check failed.'
        receipt = dict(status='PASS' if error is None else 'FAIL', task=args.task,
                       mode=args.mode if args.task == 'quadrature' else None,
                       case=args.case if args.task == 'packet' else None,
                       commands=commands, task_result=task_result,
                       frozen_repository_unchanged=post.returncode == 0,
                       wall_seconds=time.perf_counter() - start, error=error,
                       environment=dict(python=sys.version, executable=sys.executable,
                                        child_optimization='0', threads=env['OMP_NUM_THREADS']))
        (out / 'REPLAY_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
    return 0 if error is None else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print('Replay error: ' + str(exc), file=sys.stderr)
        sys.exit(1)
