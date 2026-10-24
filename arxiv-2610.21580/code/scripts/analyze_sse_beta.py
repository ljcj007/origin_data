#!/usr/bin/env python3
"""汇总 SSE 的 beta/L^2 收敛扫描。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def load_rows(run_dir: Path) -> list[dict[str, float | int | str]]:
    spec = json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))
    expected_source = spec["provenance"]["source_sha256"].upper()
    expected_executable = spec["provenance"]["executable_sha256"].upper()
    rows: list[dict[str, float | int | str]] = []
    for cell in spec["cells"]:
        path = run_dir / "cells" / cell["cell_id"] / "manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        provenance = manifest["provenance"]
        if not manifest["status"]["success"]:
            raise RuntimeError(f"失败单元：{cell['cell_id']}")
        if provenance["source_sha256"].upper() != expected_source:
            raise RuntimeError(f"源码哈希不一致：{cell['cell_id']}")
        if provenance["executable_sha256"].upper() != expected_executable:
            raise RuntimeError(f"程序哈希不一致：{cell['cell_id']}")

        params = manifest["params"]
        obs = manifest["observables"]
        beta = float(params["beta"])
        winding = float(obs["winding_squared"]["value"])
        winding_error = float(obs["winding_squared"]["error"])
        target = float(manifest["benchmarks"]["zero_temperature_ferrimagnetic_m2"])
        ferri = float(obs["ferrimagnetic_m2"]["value"])
        ferri_error = float(obs["ferrimagnetic_m2"]["error"])
        rows.append(
            {
                "cell_id": cell["cell_id"],
                "L": int(params["L"]),
                "beta_over_L2": float(params["beta_over_L2"]),
                "beta": beta,
                "energy": float(obs["energy_per_site"]["value"]),
                "energy_error": float(obs["energy_per_site"]["error"]),
                "ferrimagnetic_m2": ferri,
                "ferrimagnetic_m2_error": ferri_error,
                "lieb_mattis_target": target,
                "ferrimagnetic_ratio": ferri / target,
                "ferrimagnetic_ratio_error": ferri_error / target,
                "winding_squared": winding,
                "winding_squared_error": winding_error,
                "stiffness_cell_average_corrected": 0.75 * winding / beta,
                "stiffness_cell_average_corrected_error": 0.75 * winding_error / beta,
                "lag1_energy": float(manifest["diagnostics"]["energy_bin_lag1"]),
                "lag1_ferrimagnetic_m2": float(
                    manifest["diagnostics"]["ferrimagnetic_m2_bin_lag1"]
                ),
            }
        )
    return sorted(rows, key=lambda row: (int(row["L"]), float(row["beta_over_L2"])))


def z_difference(a: dict, b: dict, key: str, error_key: str) -> float:
    sigma = math.hypot(float(a[error_key]), float(b[error_key]))
    return abs(float(a[key]) - float(b[key])) / sigma if sigma else math.inf


def convergence_summary(rows: list[dict[str, float | int | str]]) -> dict:
    sizes = sorted({int(row["L"]) for row in rows})
    checks = []
    for size in sizes:
        series = [row for row in rows if int(row["L"]) == size]
        previous, final = series[-2], series[-1]
        target_sigma = abs(
            float(final["ferrimagnetic_m2"]) - float(final["lieb_mattis_target"])
        ) / float(final["ferrimagnetic_m2_error"])
        checks.append(
            {
                "L": size,
                "largest_beta_over_L2": float(final["beta_over_L2"]),
                "energy_change_sigma_from_previous": z_difference(
                    previous, final, "energy", "energy_error"
                ),
                "ferrimagnetic_change_sigma_from_previous": z_difference(
                    previous,
                    final,
                    "ferrimagnetic_m2",
                    "ferrimagnetic_m2_error",
                ),
                "ferrimagnetic_distance_from_Lieb_Mattis_sigma": target_sigma,
            }
        )
    return {
        "cell_count": len(rows),
        "max_abs_energy_bin_lag1": max(abs(float(row["lag1_energy"])) for row in rows),
        "max_abs_ferrimagnetic_m2_bin_lag1": max(
            abs(float(row["lag1_ferrimagnetic_m2"])) for row in rows
        ),
        "largest_ratio_checks": checks,
        "recommendation": {
            "L_ge_4": "beta=2 L^2 可作为第一版基态投影长度",
            "L_2": "beta/L^2=2 尚未收敛，建议至少补到 4",
            "stiffness": "由原始绕数按 0.75<W^2>/beta 后处理修正",
        },
    }


def write_csv(rows: list[dict[str, float | int | str]], output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_svg(rows: list[dict[str, float | int | str]], output: Path) -> None:
    width, height = 920, 760
    left, right, top, gap, bottom = 100, 42, 42, 70, 76
    panel_height = (height - top - bottom - gap) / 2
    plot_width = width - left - right
    ratios = [float(row["beta_over_L2"]) for row in rows]
    energies = [float(row["energy"]) for row in rows]
    ferri_ratios = [float(row["ferrimagnetic_ratio"]) for row in rows]
    xmin, xmax = min(ratios) - 0.05, max(ratios) + 0.05
    energy_min, energy_max = min(energies) - 0.012, max(energies) + 0.012
    ferri_min, ferri_max = min(ferri_ratios) - 0.08, max(ferri_ratios) + 0.08
    sizes = sorted({int(row["L"]) for row in rows})
    colors = ["#2166ac", "#1b9e77", "#d95f02", "#7570b3"]

    def x(value: float) -> float:
        return left + (value - xmin) * plot_width / (xmax - xmin)

    def y(value: float, low: float, high: float, panel_top: float) -> float:
        return panel_top + (high - value) * panel_height / (high - low)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.tick{font-size:14px}.label{font-size:18px}.legend{font-size:14px}.panel{font-size:19px;font-weight:bold}</style>',
    ]
    panel_tops = [top, top + panel_height + gap]
    for panel_top in panel_tops:
        parts.append(
            f'<rect x="{left}" y="{panel_top:.2f}" width="{plot_width}" height="{panel_height:.2f}" fill="none" stroke="#222" stroke-width="1.4"/>'
        )
        for ratio in sorted(set(ratios)):
            xx = x(ratio)
            parts.append(
                f'<line x1="{xx:.2f}" y1="{panel_top:.2f}" x2="{xx:.2f}" y2="{panel_top+panel_height:.2f}" stroke="#eeeeee"/>'
            )

    for panel_top, low, high, digits in [
        (panel_tops[0], energy_min, energy_max, 3),
        (panel_tops[1], ferri_min, ferri_max, 2),
    ]:
        for index in range(5):
            value = low + index * (high - low) / 4
            yy = y(value, low, high, panel_top)
            parts.append(
                f'<line x1="{left}" y1="{yy:.2f}" x2="{width-right}" y2="{yy:.2f}" stroke="#eeeeee"/>'
            )
            parts.append(
                f'<text class="tick" x="{left-12}" y="{yy+5:.2f}" text-anchor="end">{value:.{digits}f}</text>'
            )

    for index, size in enumerate(sizes):
        color = colors[index % len(colors)]
        series = [row for row in rows if int(row["L"]) == size]
        for panel, (key, error_key, low, high) in enumerate(
            [
                ("energy", "energy_error", energy_min, energy_max),
                ("ferrimagnetic_ratio", "ferrimagnetic_ratio_error", ferri_min, ferri_max),
            ]
        ):
            panel_top = panel_tops[panel]
            points = " ".join(
                f'{x(float(row["beta_over_L2"])):.2f},{y(float(row[key]), low, high, panel_top):.2f}'
                for row in series
            )
            parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>')
            for row in series:
                xx = x(float(row["beta_over_L2"]))
                yy = y(float(row[key]), low, high, panel_top)
                y1 = y(float(row[key]) + float(row[error_key]), low, high, panel_top)
                y2 = y(float(row[key]) - float(row[error_key]), low, high, panel_top)
                parts.append(f'<line x1="{xx:.2f}" y1="{y1:.2f}" x2="{xx:.2f}" y2="{y2:.2f}" stroke="{color}"/>')
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.8" fill="white" stroke="{color}" stroke-width="2"/>')
        legend_x, legend_y = width - right - 95, top + 25 + 23 * index
        parts.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x+22}" y2="{legend_y}" stroke="{color}" stroke-width="2"/>')
        parts.append(f'<text class="legend" x="{legend_x+29}" y="{legend_y+5}">L={size}</text>')

    target_y = y(1.0, ferri_min, ferri_max, panel_tops[1])
    parts.append(
        f'<line x1="{left}" y1="{target_y:.2f}" x2="{width-right}" y2="{target_y:.2f}" stroke="#555" stroke-dasharray="7 5"/>'
    )
    parts.append(f'<text class="panel" x="{left+10}" y="{top+26}">(a)</text>')
    parts.append(f'<text class="panel" x="{left+10}" y="{panel_tops[1]+26}">(b)</text>')
    parts.append(f'<text class="label" x="{left+plot_width/2:.2f}" y="{height-22}" text-anchor="middle">β / L²</text>')
    parts.append(f'<text class="label" transform="translate(31 {top+panel_height/2:.2f}) rotate(-90)" text-anchor="middle">Energy per site</text>')
    parts.append(f'<text class="label" transform="translate(31 {panel_tops[1]+panel_height/2:.2f}) rotate(-90)" text-anchor="middle">m_f^2 / Lieb-Mattis target</text>')
    for ratio in sorted(set(ratios)):
        parts.append(f'<text class="tick" x="{x(ratio):.2f}" y="{panel_tops[1]+panel_height+27:.2f}" text-anchor="middle">{ratio:g}</text>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = load_rows(args.run_dir)
    summary = convergence_summary(rows)
    write_csv(rows, args.run_dir / "scan_observables.csv")
    (args.run_dir / "beta_convergence_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_svg(rows, args.run_dir / "beta_convergence.svg")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
