"""Spot-check the report's bandwidth screen with the supplied independent solver.
Compare all 193 settings at frequencies with both saved Green blocks and comparison rows.
"""
import sys, json
from pathlib import Path
import numpy as np
from alpha_solver import solve_alpha, load_g0

fs = json.load(open(Path(__file__).resolve().parents[1] / 'comparison_inputs' / 'FINAL_SCREEN.json'))
surf = {(r['omega'], r['setting']): r for r in fs['surface_rows']}
tet = {(r['omega'], r['setting']): r for r in fs['independent_rows']}
omegas = [float(x) for x in sys.argv[1:]] or [1.1125, 1.3]
for om in omegas:
    comparison_rows = [r for (w, s), r in surf.items() if abs(w - om) < 1e-9]
    if not comparison_rows:
        print(f"omega={om}: NO COMPARISON DATA in FINAL_SCREEN.json; no validation claim made")
        continue
    for Np, deg in ((2048, 2), (4096, 2)):
        try:
            g0, ImT, ImL = load_g0(om, Np=Np, deg=deg)
        except KeyError:
            continue
        rows = comparison_rows
        mine = []
        for r in rows:
            out = solve_alpha(np.array(r['g']), om, g0, ImT, ImL, 1.69, alpha=1.0)
            t = tet.get((r['omega'], r['setting']))
            mine.append((r['setting'], out['TqT'], r['TqT'], t['TT_modal'] if t else np.nan, out['closure']))
        d_s = np.array([m[1] - m[2] for m in mine]); d_t = np.array([m[1] - m[3] for m in mine])
        imin = int(np.argmin([m[1] for m in mine]))
        print(f'omega={om} Re N={Np}: {len(mine)} settings; mine - surface: mean {d_s.mean():+.2e} max|.| {np.abs(d_s).max():.2e}; '
              f'mine - tetra: mean {np.nanmean(d_t):+.2e} max|.| {np.nanmax(np.abs(d_t)):.2e}; max|closure| {max(abs(m[4]) for m in mine):.1e}')
        print(f'   my minimum: {mine[imin][0]} TqT={mine[imin][1]:.10f} (surface {mine[imin][2]:.10f}, tetra {mine[imin][3]:.10f}); '
              f'nominal: {[m for m in mine if m[0] == "nominal"][0][1]:.10f}')
