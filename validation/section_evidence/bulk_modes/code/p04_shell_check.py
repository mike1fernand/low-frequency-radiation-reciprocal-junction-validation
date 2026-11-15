from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

from p04_common import (
    TARGET,
    eig,
    group_velocity,
    jsonable,
    qt_projector,
    shell_point,
    source_full,
    source_radial,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", type=int, choices=[0, 1, 2], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    branch = args.branch

    def objective(angles: np.ndarray, kind: str) -> float:
        k, vals, vecs = shell_point(branch, angles)
        khat = k / np.linalg.norm(k)
        p_qt = qt_projector(vecs)
        if kind == "cluster_leak":
            return -float(khat @ p_qt @ khat)
        if kind == "own":
            own = float(abs(khat @ vecs[:, branch]) ** 2)
            return -own if branch < 2 else own
        if kind == "gap":
            return float(vals[2] - vals[1])
        raise ValueError(kind)

    kinds = ["cluster_leak", "own"] if branch < 2 else ["own", "gap"]
    results = {}
    for index, kind in enumerate(kinds):
        result = differential_evolution(
            lambda x: objective(x, kind),
            [(0.0, math.pi), (0.0, 2.0 * math.pi)],
            seed=20260813 + 10 * branch + index,
            tol=1e-10,
            maxiter=180,
            popsize=16,
            polish=True,
            workers=1,
        )
        k, vals, vecs = shell_point(branch, result.x)
        khat = k / np.linalg.norm(k)
        p_qt = qt_projector(vecs)
        results[kind] = {
            "success": bool(result.success),
            "optimizer_message": str(result.message),
            "raw_objective": float(result.fun),
            "k": k,
            "shell_residual": float(vals[branch] - TARGET),
            "cluster_leak": float(khat @ p_qt @ khat),
            "own": float(abs(khat @ vecs[:, branch]) ** 2),
            "gap": float(vals[2] - vals[1]),
            "group_speed": float(np.linalg.norm(group_velocity(k, branch))),
            "nfev": int(result.nfev),
        }

    witness_points = {
        0: np.array([1.022121111600022, 0.7119903394245296, -0.020431525569168604]),
        1: np.array([0.8807183236959493, 0.21703145119419856, -0.8514466270114215]),
    }
    witness = None
    if branch < 2:
        k = witness_points[branch]
        vals, vecs = eig(k)
        f_radial = source_radial(k)
        f_full = source_full(k)
        radial_sq = float(abs(np.vdot(vecs[:, branch], f_radial)) ** 2)
        full_sq = float(abs(np.vdot(vecs[:, branch], f_full)) ** 2)
        witness = {
            "source_definition": "complete bare finite-slab X+R line-driven source; radial component also reported",
            "k": k,
            "shell_residual": float(vals[branch] - TARGET),
            "group_speed": float(np.linalg.norm(group_velocity(k, branch))),
            "lower_gap": float(abs(vals[1] - vals[0])),
            "upper_gap": float(vals[2] - vals[1]),
            "radial_overlap_sq": radial_sq,
            "full_overlap_sq": full_sq,
            "radial_normalized_overlap": radial_sq / max(float(np.vdot(f_radial, f_radial).real), 1e-30),
            "full_normalized_overlap": full_sq / max(float(np.vdot(f_full, f_full).real), 1e-30),
        }

    checks = {
        "optimization_success": all(item["success"] for item in results.values()),
        "shell_residual": max(abs(item["shell_residual"]) for item in results.values()) < 1e-12,
    }
    if branch < 2:
        checks["cluster_bound"] = results["cluster_leak"]["cluster_leak"] < 1.25e-3
        checks["regular_positive_witness"] = bool(
            witness is not None
            and abs(witness["shell_residual"]) < 5e-12
            and witness["group_speed"] > 0.1
            and witness["lower_gap"] > 1e-3
            and witness["upper_gap"] > 0.3
            and witness["full_overlap_sq"] > 1e-2
        )
    else:
        checks["purity_and_gap"] = bool(
            results["own"]["own"] > 0.99988 and results["gap"]["gap"] > 0.9390
        )

    status = "passed" if all(checks.values()) else "failed"
    output = {
        "schema": "paper1.p04.shell-check.v1",
        "status": status,
        "branch": branch,
        "checks": checks,
        "results": results,
        "witness": witness,
        "producer_contract": {
            "fixed_random_seeds": [20260813 + 10 * branch, 20260814 + 10 * branch],
            "deterministic_inputs": True,
            "contains_timestamps": False,
            "external_reproduction_claim": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, default=jsonable) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "branch": branch, "checks": checks}, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
