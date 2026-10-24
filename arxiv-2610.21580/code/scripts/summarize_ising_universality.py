#!/usr/bin/env python3
"""汇总经典 Ising 普适类一致性拟合。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIT_SHA256 = "9D6C8ADED7A388AE15AFE84C5B5BF2C3324573162E6C62E731E4B2FE4039FB32"


def read_fit(path: Path) -> dict[str, dict[str, float | None]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result: dict[str, dict[str, float | None]] = {}
    for row in rows:
        result[row["parameter"]] = {
            key: float(row[key]) if row[key] else None for key in ("value", "ci_lo", "ci_hi")
        }
    if int(result["failed_resamples"]["value"] or 0) != 0:
        raise RuntimeError(f"bootstrap 存在失败：{path}")
    return result


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def power_summary(fit: dict[str, dict[str, float | None]], growth: bool) -> dict[str, float | bool]:
    alpha = fit["alpha"]
    sign = -1.0 if growth else 1.0
    value = sign * float(alpha["value"])
    ci_values = [sign * float(alpha["ci_lo"]), sign * float(alpha["ci_hi"])]
    expected = 1.75 if growth else 0.125
    return {
        "value": value,
        "ci95_low": min(ci_values),
        "ci95_high": max(ci_values),
        "expected_2d_ising": expected,
        "expected_inside_statistical_ci": min(ci_values) <= expected <= max(ci_values),
        "chi2_per_dof": float(fit["chi2_per_dof"]["value"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir", type=Path, default=ROOT / "results" / "ising_tc_wolff_scaling_v1"
    )
    parser.add_argument(
        "--fit-script",
        type=Path,
        default=Path(
            "D:/Project/harnessing-quantum-fork/submission/quantum.harness/scripts/scaling_fit.py"
        ),
    )
    args = parser.parse_args()
    if sha256(args.fit_script) != FIT_SHA256:
        raise RuntimeError("拟合程序哈希不匹配")

    magnetization = read_fit(args.run_dir / "abs_m_power.csv")
    susceptibility = read_fit(args.run_dir / "connected_chi_power.csv")
    magnetization_lmin16 = read_fit(args.run_dir / "abs_m_power_Lmin16.csv")
    susceptibility_lmin16 = read_fit(args.run_dir / "connected_chi_power_Lmin16.csv")
    binder = read_fit(args.run_dir / "binder_collapse.csv")
    summary = {
        "schema_version": 1,
        "status": "consistency_check_only",
        "fit_script_sha256": FIT_SHA256,
        "fixed_exact_temperature": 2.0 / math.acosh(2.0 + math.sqrt(3.0)),
        "primary_Lmin8": {
            "abs_staggered_m_decay_exponent": power_summary(magnetization, growth=False),
            "connected_staggered_susceptibility_growth_exponent": power_summary(
                susceptibility, growth=True
            ),
        },
        "sensitivity_Lmin16": {
            "abs_staggered_m_decay_exponent": power_summary(magnetization_lmin16, growth=False),
            "connected_staggered_susceptibility_growth_exponent": power_summary(
                susceptibility_lmin16, growth=True
            ),
        },
        "binder_collapse": {
            "nu_fit_label": float(binder["nu"]["value"]),
            "ci95_low": float(binder["nu"]["ci_lo"]),
            "ci95_high": float(binder["nu"]["ci_hi"]),
            "collapse_residual": float(binder["collapse_residual"]["value"]),
            "interpretation": "描述性检查；分箱目标函数有平台，不作为独立精度测量",
        },
        "cross_method_check": {
            "status": "statistical_agreement_but_full_budget_missing",
            "reason": "精确二维 Ising 指数落入统计区间，但温度插值、有限尺寸修正和单热化历史没有完整数值预算",
            "allowed_claim": "Wolff 数据与二维 Ising 普适类一致",
            "forbidden_claim": "由该扫描独立精密测定临界指数",
        },
    }
    output = args.run_dir / "fit_summary.json"
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"写出 {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
