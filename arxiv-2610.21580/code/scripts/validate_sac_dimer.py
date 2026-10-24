#!/usr/bin/env python3
"""用独立 Heisenberg 二聚体检查虚时关联与 SAC 峰位。"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def read_tq(path: Path, mode: int) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()
    rows: list[list[float]] = []
    active = False
    for line in lines:
        fields = line.split()
        if len(fields) == 1:
            active = int(fields[0]) == mode
        elif active and len(fields) == 3:
            rows.append([float(value) for value in fields])
    if not rows:
        raise ValueError(f"tq.dat 中没有 mode={mode}")
    return np.asarray(rows)


def exact_normalized_correlation(beta: float, tau: np.ndarray, exchange: float = 1.0) -> np.ndarray:
    numerator = np.exp(-exchange * tau) + np.exp(-exchange * (beta - tau))
    return numerator / (1.0 + math.exp(-beta * exchange))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--mode", type=int, default=2)
    parser.add_argument("--beta", type=float, default=8.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    tq = read_tq(run_dir / "tq.dat", args.mode)
    normalized = tq[:, 1] / tq[0, 1]
    exact = exact_normalized_correlation(args.beta, tq[:, 0])
    imaginary_time_max_error = float(np.max(np.abs(normalized - exact)))

    spectrum = np.loadtxt(run_dir / "sw.dat")
    omega = spectrum[:, 0]
    weight = spectrum[:, 1]
    peak = float(omega[int(np.argmax(weight))])
    integral = float(np.trapz(weight, omega))
    mean = float(np.trapz(omega * weight, omega) / integral)
    expected_c0 = 1.5 * (math.exp(args.beta) + 1.0) / (math.exp(args.beta) + 3.0)
    expected_integral = math.pi * expected_c0

    checks = {
        "imaginary_time_max_abs_error": imaginary_time_max_error < 0.03,
        "sac_peak": abs(peak - 1.0) < 0.15,
        "sac_first_moment": abs(mean - 1.0) < 0.15,
        "spectral_sum": abs(integral - expected_integral) / expected_integral < 0.05,
    }
    payload = {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "beta": args.beta,
        "exact_gap": 1.0,
        "imaginary_time_max_abs_error": imaginary_time_max_error,
        "sac_peak": peak,
        "sac_first_moment": mean,
        "spectral_integral": integral,
        "expected_spectral_integral": expected_integral,
    }
    output = args.output or run_dir / "validation.json"
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    raise SystemExit(0 if payload["status"] == "passed" else 1)


if __name__ == "__main__":
    main()
