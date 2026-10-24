#!/usr/bin/env python3
"""整理 SSE 外推、截断敏感性和 Lieb--Mattis 检验。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


ENERGY_COEFFICIENTS = ("c0", "c1", "c2", "c3", "c4", "c5", "c6")
ENERGY_PINS = ("c1", "c2", "c3", "c5")
EVEN_COEFFICIENTS = ("c0", "c1", "c2", "c3", "c4")
EVEN_PINS = ("c1", "c3")
EXPECTED_SIZES = (4, 6, 8, 10, 12, 14, 16)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"空数据表：{path}")
    return rows


def finite_float(value: str, label: str, path: Path) -> float:
    """读取有限浮点数。"""
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} 不是数值：{path}") from exc
    if not math.isfinite(result):
        raise RuntimeError(f"{label} 不是有限值：{path}")
    return result


def read_fit(
    path: Path,
    required: tuple[str, ...],
    pinned_zero: tuple[str, ...],
) -> dict[str, Any]:
    rows = read_rows(path)
    parameters = [row.get("parameter", "") for row in rows]
    if len(parameters) != len(set(parameters)):
        raise RuntimeError(f"拟合输出存在重复参数：{path}")
    values = {row["parameter"]: row for row in rows}
    names = (*required, "chi2_per_dof", "failed_resamples")
    missing = [name for name in names if name not in values]
    if missing:
        raise RuntimeError(f"拟合输出不完整：{path}")
    coefficients = {name: finite_float(values[name]["value"], name, path) for name in required}
    coefficient_ci95: dict[str, list[float]] = {}
    for name in required:
        lower = finite_float(values[name].get("ci_lo", ""), f"{name} ci_lo", path)
        upper = finite_float(values[name].get("ci_hi", ""), f"{name} ci_hi", path)
        if lower > coefficients[name] or coefficients[name] > upper:
            raise RuntimeError(f"{name} 不在 95% CI 内：{path}")
        coefficient_ci95[name] = [lower, upper]
    intercept = coefficients["c0"]
    failed_value = finite_float(values["failed_resamples"]["value"], "failed_resamples", path)
    if not failed_value.is_integer() or failed_value < 0:
        raise RuntimeError(f"failed_resamples 不是非负整数：{path}")
    failed = int(failed_value)
    if failed:
        raise RuntimeError(f"bootstrap 存在 {failed} 个失败样本：{path}")
    invalid_pins = [name for name in pinned_zero if abs(coefficients[name]) > 1.0e-12]
    if invalid_pins:
        raise RuntimeError(f"预注册固定系数不是零 {invalid_pins}：{path}")
    chi2 = finite_float(values["chi2_per_dof"]["value"], "chi2_per_dof", path)
    if chi2 < 0.0:
        raise RuntimeError(f"chi2_per_dof 不能为负：{path}")
    return {
        "intercept": intercept,
        "ci95": coefficient_ci95["c0"],
        "coefficients": coefficients,
        "coefficient_ci95": coefficient_ci95,
        "chi2_per_dof": chi2,
        "failed_resamples": failed,
        "source": str(path),
    }


def validate_fss_table(rows: list[dict[str, str]]) -> None:
    """锁定七尺寸正式表及拟合误差。"""
    required = (
        "L",
        "N",
        "beta",
        "replicas",
        "pooled_bins",
        "energy_per_site",
        "energy_per_site_error",
        "staggered_m2",
        "staggered_m2_error",
        "ferrimagnetic_m2",
        "ferrimagnetic_m2_error",
        "lieb_mattis_target",
        "sublattice_plus_m2",
        "sublattice_plus_m2_error",
        "sublattice_minus_m2",
        "sublattice_minus_m2_error",
        "sublattice_cross",
        "sublattice_cross_error",
        "spin_stiffness_cell_average",
        "spin_stiffness_cell_average_error",
    )
    missing = [name for name in required if name not in rows[0]]
    if missing:
        raise RuntimeError(f"FSS 表缺少列：{missing}")
    if len(rows) != len(EXPECTED_SIZES):
        raise RuntimeError("FSS 表必须恰含七个尺寸")
    sizes: list[int] = []
    for row in rows:
        numeric = {name: finite_float(row[name], name, Path("fss_observables.csv")) for name in required}
        size = int(numeric["L"])
        if numeric["L"] != size:
            raise RuntimeError("FSS 尺寸必须为整数")
        sizes.append(size)
        expected = {
            "N": 5 * size * size,
            "beta": 2 * size * size,
            "replicas": 2,
            "pooled_bins": 128,
        }
        for name, target in expected.items():
            if numeric[name] != target:
                raise RuntimeError(f"FSS 协议字段不一致：L={size}, {name}")
        for name in required:
            if name.endswith("_error") and numeric[name] <= 0.0:
                raise RuntimeError(f"FSS 误差必须为正：L={size}, {name}")
    if sizes != list(EXPECTED_SIZES):
        raise RuntimeError("FSS 尺寸轴不完整或未排序")


def validate_fit_against_table(
    fit: dict[str, Any],
    rows: list[dict[str, str]],
    observable: str,
    coefficients: tuple[str, ...],
    pinned_zero: tuple[str, ...],
) -> None:
    """重算残差，确认拟合使用了声明列和系数。"""
    error_name = f"{observable}_error"
    if observable not in rows[0] or error_name not in rows[0]:
        raise RuntimeError(f"FSS 表缺少拟合列：{observable}")
    chi2 = 0.0
    for row in rows:
        size = finite_float(row["L"], "L", Path("fss_observables.csv"))
        value = finite_float(row[observable], observable, Path("fss_observables.csv"))
        error = finite_float(row[error_name], error_name, Path("fss_observables.csv"))
        if error <= 0.0:
            raise RuntimeError(f"拟合误差必须为正：{observable}, L={size:g}")
        prediction = sum(
            float(fit["coefficients"][name]) / size ** int(name[1:]) for name in coefficients
        )
        chi2 += ((value - prediction) / error) ** 2
    free_count = len(coefficients) - len(pinned_zero)
    expected = chi2 / max(len(rows) - free_count, 1)
    if not math.isclose(expected, float(fit["chi2_per_dof"]), rel_tol=1.0e-9, abs_tol=1.0e-10):
        raise RuntimeError(f"拟合输出与声明数据列不一致：{observable}")


def validate_diagnostics(diagnostics: dict[str, Any]) -> None:
    """核对汇总文件的正式运行规模。"""
    expected = {
        "cell_count": 14,
        "size_count": 7,
        "sizes": list(EXPECTED_SIZES),
        "replicas_per_size": 2,
        "bins_per_size": [128] * 7,
    }
    for name, target in expected.items():
        if diagnostics.get(name) != target:
            raise RuntimeError(f"FSS 汇总协议字段不一致：{name}")


def transform_sublattice(staggered_m2: float) -> dict[str, float]:
    if not math.isfinite(staggered_m2) or staggered_m2 <= 0.0:
        raise RuntimeError("staggered_m2 外推值必须为正的有限数")
    staggered = math.sqrt(staggered_m2)
    plus = (5.0 * staggered - 0.5) / 4.0
    minus = (5.0 * staggered + 0.5) / 6.0
    return {
        "staggered_moment": staggered,
        "plus_moment": plus,
        "minus_moment": minus,
        "plus_m2": plus * plus,
        "minus_m2": minus * minus,
        "cross": -plus * minus,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    fits = run_dir / "fits"

    primary_files = {
        "energy_per_site": ("energy_fit.csv", ENERGY_COEFFICIENTS, ENERGY_PINS),
        "staggered_m2": ("staggered_m2_fit.csv", EVEN_COEFFICIENTS, EVEN_PINS),
        "spin_stiffness_cell_average": ("stiffness_fit.csv", EVEN_COEFFICIENTS, EVEN_PINS),
    }
    secondary_files = {
        name: f"{name}_secondary_fit.csv"
        for name in ("sublattice_plus_m2", "sublattice_minus_m2", "sublattice_cross")
    }
    primary = {
        name: read_fit(fits / filename, required, pins)
        for name, (filename, required, pins) in primary_files.items()
    }
    secondary = {
        name: read_fit(fits / filename, EVEN_COEFFICIENTS, EVEN_PINS)
        for name, filename in secondary_files.items()
    }

    sensitivity: dict[str, dict[str, Any]] = {}
    sensitivity_names = {
        "energy_per_site": "energy",
        "staggered_m2": "staggered_m2",
        "spin_stiffness_cell_average": "spin_stiffness_cell_average",
    }
    for observable, stem in sensitivity_names.items():
        sensitivity[observable] = {
            f"Lmin_{minimum}": read_fit(
                fits / "sensitivity" / f"{stem}_Lmin{minimum}.csv",
                ENERGY_COEFFICIENTS if observable == "energy_per_site" else EVEN_COEFFICIENTS,
                ENERGY_PINS if observable == "energy_per_site" else EVEN_PINS,
            )
            for minimum in (6, 8)
        }
        intercepts = [primary[observable]["intercept"]] + [
            item["intercept"] for item in sensitivity[observable].values()
        ]
        sensitivity[observable]["intercept_span"] = [min(intercepts), max(intercepts)]

    table = read_rows(run_dir / "fss_observables.csv")
    validate_fss_table(table)
    for observable, (_, coefficients, pins) in primary_files.items():
        validate_fit_against_table(primary[observable], table, observable, coefficients, pins)
    for observable in secondary_files:
        validate_fit_against_table(
            secondary[observable], table, observable, EVEN_COEFFICIENTS, EVEN_PINS
        )
    for observable, stem in sensitivity_names.items():
        coefficients = ENERGY_COEFFICIENTS if observable == "energy_per_site" else EVEN_COEFFICIENTS
        pins = ENERGY_PINS if observable == "energy_per_site" else EVEN_PINS
        for minimum in (6, 8):
            subset = [row for row in table if int(row["L"]) >= minimum]
            validate_fit_against_table(
                sensitivity[observable][f"Lmin_{minimum}"],
                subset,
                observable,
                coefficients,
                pins,
            )
    ferrimagnetic_residuals = []
    for row in table:
        value = finite_float(row["ferrimagnetic_m2"], "ferrimagnetic_m2", run_dir)
        target = finite_float(row["lieb_mattis_target"], "lieb_mattis_target", run_dir)
        error = finite_float(row["ferrimagnetic_m2_error"], "ferrimagnetic_m2_error", run_dir)
        residual = (value - target) / error
        ferrimagnetic_residuals.append({"L": int(row["L"]), "standardized_residual": residual})

    staggered_fit = primary["staggered_m2"]
    constrained = transform_sublattice(staggered_fit["intercept"])
    constrained_lo = transform_sublattice(staggered_fit["ci95"][0])
    constrained_hi = transform_sublattice(staggered_fit["ci95"][1])
    constrained["ci95"] = {
        key: sorted([constrained_lo[key], constrained_hi[key]])
        for key in constrained_lo
    }

    diagnostics = json.loads((run_dir / "fss_summary.json").read_text(encoding="utf-8"))
    validate_diagnostics(diagnostics)
    result = {
        "uncertainty_boundary": (
            "输入误差取未分块及 replica-aware 2/4/8-bin 分块标准误的最大值；区间来自固定拟合形式下的参数 "
            "bootstrap，不包含拟合形式或 beta*J/L^2 截断的总系统误差"
        ),
        "primary": primary,
        "secondary": secondary,
        "constrained_sublattice_assuming_collinear_order_and_cross_factorization": constrained,
        "sensitivity": sensitivity,
        "ferrimagnetic_exact_check": {
            "target": "0.01+0.02/L^2",
            "residuals": ferrimagnetic_residuals,
            "max_abs_standardized_residual": max(
                abs(item["standardized_residual"]) for item in ferrimagnetic_residuals
            ),
        },
        "bin_diagnostics": {
            "max_abs_replica_lag1": diagnostics["max_abs_replica_lag1"],
            "max_blocking_error_ratio": diagnostics["max_blocking_error_ratio"],
            "max_replica_difference_sigma": diagnostics["max_replica_difference_sigma"],
        },
    }
    output = run_dir / "fit_results_summary.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
