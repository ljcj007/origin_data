#!/usr/bin/env python3
"""从 Wolff 原始 bin 构造临界点有限尺寸观测量。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
EXACT_TC = 2.0 / math.acosh(2.0 + math.sqrt(3.0))
BLOCK_SIZES = (1, 2, 4, 8)


def read_bins(path: Path) -> dict[str, list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 64:
        raise RuntimeError(f"bin 数量不是 64：{path}")
    return {name: [float(row[name]) for row in rows] for name in rows[0]}


def jackknife_error(
    columns: dict[str, list[float]],
    block_size: int,
    estimator: Callable[[dict[str, list[float]]], float],
) -> float:
    count = len(next(iter(columns.values())))
    if count % block_size:
        raise ValueError("bin 数量不能被分块大小整除")
    replicas: list[float] = []
    for start in range(0, count, block_size):
        stop = start + block_size
        reduced = {name: values[:start] + values[stop:] for name, values in columns.items()}
        replicas.append(estimator(reduced))
    center = statistics.fmean(replicas)
    blocks = len(replicas)
    return math.sqrt((blocks - 1.0) / blocks * sum((value - center) ** 2 for value in replicas))


def summarize_cell(cell_dir: Path) -> dict[str, float | int]:
    manifest = json.loads((cell_dir / "manifest.json").read_text(encoding="utf-8"))
    if not manifest.get("status", {}).get("success"):
        raise RuntimeError(f"失败单元：{cell_dir.name}")
    size = int(manifest["params"]["L"])
    temperature = float(manifest["params"]["temperature"])
    bins = read_bins(cell_dir / "ising_bins.csv")
    columns = {
        "abs_m": bins["abs_staggered_m"],
        "m2": bins["staggered_m2"],
    }
    nsite = 5 * size * size

    def abs_m(data: dict[str, list[float]]) -> float:
        return statistics.fmean(data["abs_m"])

    def connected_chi(data: dict[str, list[float]]) -> float:
        return nsite / temperature * (
            statistics.fmean(data["m2"]) - statistics.fmean(data["abs_m"]) ** 2
        )

    abs_value = abs_m(columns)
    chi_value = connected_chi(columns)
    abs_errors = [jackknife_error(columns, block, abs_m) for block in BLOCK_SIZES]
    chi_errors = [jackknife_error(columns, block, connected_chi) for block in BLOCK_SIZES]
    recorded_chi = statistics.fmean(bins["connected_staggered_susceptibility"])
    return {
        "L": size,
        "temperature": temperature,
        "abs_m": abs_value,
        "abs_m_error": max(abs_errors),
        "connected_chi": chi_value,
        "connected_chi_error": max(chi_errors),
        "recorded_binwise_connected_chi": recorded_chi,
        "binwise_bias": recorded_chi - chi_value,
    }


def interpolate(low: dict[str, float | int], high: dict[str, float | int]) -> dict[str, float | int]:
    low_t = float(low["temperature"])
    high_t = float(high["temperature"])
    weight = (EXACT_TC - low_t) / (high_t - low_t)
    row: dict[str, float | int] = {"L": int(low["L"])}
    for name in ("abs_m", "connected_chi"):
        low_value = float(low[name])
        high_value = float(high[name])
        low_error = float(low[f"{name}_error"])
        high_error = float(high[f"{name}_error"])
        row[name] = (1.0 - weight) * low_value + weight * high_value
        row[f"{name}_error"] = math.hypot((1.0 - weight) * low_error, weight * high_error)
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "results" / "ising_tc_wolff_v1")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results" / "ising_tc_wolff_scaling_v1"
    )
    args = parser.parse_args()

    cells = [summarize_cell(path) for path in sorted((args.run_dir / "cells").glob("cell-*"))]
    if len(cells) != 35:
        raise RuntimeError(f"正式 Wolff 单元应为 35，实际为 {len(cells)}")
    by_size: dict[int, dict[float, dict[str, float | int]]] = {}
    for cell in cells:
        by_size.setdefault(int(cell["L"]), {})[float(cell["temperature"])] = cell
    if sorted(by_size) != [8, 12, 16, 24, 32, 48, 64]:
        raise RuntimeError("尺寸集合不符合正式协议")

    rows = [interpolate(by_size[size][1.004], by_size[size][1.008]) for size in sorted(by_size)]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "critical_observables.csv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    maximum_bias = max(cells, key=lambda row: abs(float(row["binwise_bias"])))
    summary = {
        "schema_version": 1,
        "status": "consistency_check_only",
        "exact_temperature": EXACT_TC,
        "interpolation_bracket": [1.004, 1.008],
        "interpolation_fraction": (EXACT_TC - 1.004) / 0.004,
        "sizes": sorted(by_size),
        "error_rule": "maximum contiguous delete-block jackknife error for 1,2,4,8 bins",
        "expected_2d_ising_exponent_ratios": {
            "abs_staggered_m_power": -0.125,
            "connected_staggered_susceptibility_power": 1.75,
        },
        "binwise_connected_susceptibility_boundary": {
            "reason": "每个 bin 内减去均值平方会产生有限 bin 偏差；标度表由总样本矩重构",
            "maximum_absolute_bias": abs(float(maximum_bias["binwise_bias"])),
            "signed_bias": float(maximum_bias["binwise_bias"]),
            "cell": {
                "L": int(maximum_bias["L"]),
                "temperature": float(maximum_bias["temperature"]),
            },
        },
        "systematics_not_included": [
            "0.004 温度网格上的线性插值误差",
            "有限尺寸修正",
            "单副本热化历史偏差",
        ],
    }
    (args.output_dir / "universality_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"写出 {table}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
