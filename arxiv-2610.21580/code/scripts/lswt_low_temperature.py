#!/usr/bin/env python3
"""验证二次 Goldstone 模的二维低温 LSWT 热力学。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


SITE_COUNT = 5


def geometric_coefficients(spin: float, exchange: float, lattice_constant: float) -> dict[str, float]:
    """返回物理面积、刚度和低能连续谱系数。"""
    cell_area = math.sqrt(3.0) * lattice_constant**2 / 2.0
    stiffness = exchange * spin * lattice_constant**2 / 4.0
    dos_per_area = 1.0 / (4.0 * math.pi * stiffness)
    dos_per_cell = cell_area * dos_per_area
    dos_per_site = dos_per_cell / SITE_COUNT
    energy_coefficient_per_area = math.pi / (24.0 * stiffness)
    energy_coefficient_per_cell = cell_area * energy_coefficient_per_area
    energy_coefficient_per_site = energy_coefficient_per_cell / SITE_COUNT
    heat_slope_per_area = math.pi / (12.0 * stiffness)
    heat_slope_per_cell = cell_area * heat_slope_per_area
    heat_slope_per_site = heat_slope_per_cell / SITE_COUNT
    return {
        "cell_area": cell_area,
        "stiffness": stiffness,
        "dos_per_area": dos_per_area,
        "dos_per_cell": dos_per_cell,
        "dos_per_site": dos_per_site,
        "energy_coefficient_per_area": energy_coefficient_per_area,
        "energy_coefficient_per_cell": energy_coefficient_per_cell,
        "energy_coefficient_per_site": energy_coefficient_per_site,
        "heat_slope_per_area": heat_slope_per_area,
        "heat_slope_per_cell": heat_slope_per_cell,
        "heat_slope_per_site": heat_slope_per_site,
    }


def acoustic_energy(kx: np.ndarray, ky: np.ndarray, spin: float, exchange: float) -> np.ndarray:
    """计算完整布里渊区内的声学支能量。"""
    g_squared = 3.0 + 2.0 * (np.cos(kx) + np.cos(ky) + np.cos(kx + ky))
    g = np.sqrt(np.maximum(g_squared, 0.0))
    root = np.sqrt(np.maximum(13.0 - 4.0 * g, 0.0))
    return 2.0 * exchange * spin * (3.0 - g) / (root + 1.0)


def thermal_energy_kernel(energy: np.ndarray, temperature: float) -> np.ndarray:
    """返回每个玻色模的热能。"""
    x = energy / temperature
    result = np.empty_like(x)
    small = x < 50.0
    result[small] = energy[small] / np.expm1(x[small])
    result[~small] = energy[~small] * np.exp(-x[~small])
    result[x == 0.0] = temperature
    return result


def heat_capacity_kernel(energy: np.ndarray, temperature: float) -> np.ndarray:
    """返回每个玻色模以 kB 为单位的比热。"""
    x = energy / temperature
    result = np.empty_like(x)
    small = x < 50.0
    half = 0.5 * x[small]
    result[small] = np.square(half / np.sinh(half))
    result[~small] = np.square(x[~small]) * np.exp(-x[~small])
    result[x == 0.0] = 1.0
    return result


def bz_thermodynamics(
    grid: int, temperatures_over_j: np.ndarray, spin: float, exchange: float
) -> list[dict[str, float]]:
    """用中点网格积分完整声学支。"""
    phases = -math.pi + (np.arange(grid, dtype=float) + 0.5) * (2.0 * math.pi / grid)
    energy = acoustic_energy(phases[:, None], phases[None, :], spin, exchange)
    rows: list[dict[str, float]] = []
    for temperature_over_j in temperatures_over_j:
        temperature = exchange * float(temperature_over_j)
        energy_per_site = float(np.mean(thermal_energy_kernel(energy, temperature))) / SITE_COUNT
        heat_per_site = float(np.mean(heat_capacity_kernel(energy, temperature))) / SITE_COUNT
        rows.append(
            {
                "temperature_over_J": float(temperature_over_j),
                "acoustic_energy_over_J_per_site": energy_per_site / exchange,
                "acoustic_heat_capacity_over_kB_per_site": heat_per_site,
            }
        )
    return rows


def continuum_integrals(order: int = 320, upper: float = 50.0) -> tuple[float, float]:
    """用高斯积分独立计算两个无量纲玻色积分。"""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = 0.5 * upper * (nodes + 1.0)
    scaled_weights = 0.5 * upper * weights
    energy_integrand = x / np.expm1(x)
    heat_integrand = np.square(0.5 * x / np.sinh(0.5 * x))
    return (
        float(np.dot(scaled_weights, energy_integrand)),
        float(np.dot(scaled_weights, heat_integrand)),
    )


def extrapolate_coefficients(rows: list[dict[str, float]]) -> dict[str, float]:
    """外推完整能带的低温系数并拟合幂律。"""
    temperatures = np.asarray([row["temperature_over_J"] for row in rows])
    energies = np.asarray([row["acoustic_energy_over_J_per_site"] for row in rows])
    heats = np.asarray([row["acoustic_heat_capacity_over_kB_per_site"] for row in rows])
    energy_intercept = float(np.polyfit(temperatures, energies / temperatures**2, 2)[-1])
    heat_intercept = float(np.polyfit(temperatures, heats / temperatures, 2)[-1])
    energy_exponent = float(np.polyfit(np.log(temperatures), np.log(energies), 1)[0])
    heat_exponent = float(np.polyfit(np.log(temperatures), np.log(heats), 1)[0])
    return {
        "energy_T2_coefficient_extrapolated": energy_intercept,
        "heat_T_slope_extrapolated": heat_intercept,
        "energy_loglog_exponent_over_sampled_window": energy_exponent,
        "heat_loglog_exponent_over_sampled_window": heat_exponent,
    }


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument("--exchange", type=float, default=1.0)
    parser.add_argument("--lattice-constant", type=float, default=1.0)
    parser.add_argument("--grid", type=int, default=1200)
    parser.add_argument(
        "--temperatures",
        type=float,
        nargs="+",
        default=[0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032],
        metavar="T_OVER_J",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/lswt_low_temperature_v1")
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    temperatures = np.asarray(sorted(args.temperatures), dtype=float)
    if (
        args.spin <= 0.0
        or args.exchange <= 0.0
        or args.lattice_constant <= 0.0
        or args.grid < 32
        or len(temperatures) < 4
        or np.any(temperatures <= 0.0)
    ):
        raise SystemExit("参数超出允许范围")

    geometry = geometric_coefficients(args.spin, args.exchange, args.lattice_constant)
    rows = bz_thermodynamics(args.grid, temperatures, args.spin, args.exchange)
    fits = extrapolate_coefficients(rows)
    energy_integral, heat_integral = continuum_integrals()

    dimensionless_dos_per_site = args.exchange * geometry["dos_per_site"]
    dimensionless_energy_coefficient = args.exchange * geometry["energy_coefficient_per_site"]
    dimensionless_heat_slope = args.exchange * geometry["heat_slope_per_site"]
    summary = {
        "format": "super-honeycomb-lswt-low-temperature-v1",
        "convention": {
            "dispersion": "omega = D |k_cart|^2",
            "D": geometry["stiffness"],
            "primitive_vector_length_a": args.lattice_constant,
            "cell_area": geometry["cell_area"],
            "sites_per_cell": SITE_COUNT,
            "J": args.exchange,
            "S": args.spin,
            "kB": 1.0,
        },
        "analytic_continuum": {
            "dos_per_area_per_energy": geometry["dos_per_area"],
            "dos_per_cell_per_energy": geometry["dos_per_cell"],
            "dos_per_site_per_energy": geometry["dos_per_site"],
            "J_times_dos_per_site": dimensionless_dos_per_site,
            "energy_per_site": "u/J = A (T/J)^2 + o(T^2)",
            "energy_T2_coefficient_A": dimensionless_energy_coefficient,
            "heat_per_site": "c/kB = B (T/J) + o(T)",
            "heat_T_slope_B": dimensionless_heat_slope,
        },
        "numerical_continuum_integrals": {
            "integral_x_over_expm1_x": energy_integral,
            "expected_pi_squared_over_6": math.pi**2 / 6.0,
            "integral_heat_kernel": heat_integral,
            "expected_pi_squared_over_3": math.pi**2 / 3.0,
        },
        "full_bz_acoustic_branch": {
            "midpoint_grid": args.grid,
            "temperature_points": len(rows),
            **fits,
            "energy_coefficient_relative_error": abs(
                fits["energy_T2_coefficient_extrapolated"] - dimensionless_energy_coefficient
            )
            / dimensionless_energy_coefficient,
            "heat_slope_relative_error": abs(
                fits["heat_T_slope_extrapolated"] - dimensionless_heat_slope
            )
            / dimensionless_heat_slope,
        },
        "scope": {
            "goldstone_branches": 1,
            "optical_branches_omitted": True,
            "finite_temperature_long_range_order": False,
            "reason": "二维短程各向同性 SU(2) 模型受 Mermin-Wagner 定理约束",
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "thermodynamics.csv", rows)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"u/J 的 T^2 系数：{fits['energy_T2_coefficient_extrapolated']:.10f}")
    print(f"c/kB 的 T 斜率：{fits['heat_T_slope_extrapolated']:.10f}")
    print(f"结果目录：{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
