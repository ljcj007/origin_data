#!/usr/bin/env python3
"""汇总多尺寸 Ising Binder 扫描并标记数值诊断。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def load_rows(run_dir: Path) -> list[dict[str, float | str]]:
    spec = json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))
    expected_hash = spec["provenance"]["executable_sha256"]
    rows = []
    for cell in spec["cells"]:
        manifest = json.loads(
            (run_dir / "cells" / cell["cell_id"] / "manifest.json").read_text(encoding="utf-8")
        )
        if not manifest["status"]["success"]:
            raise RuntimeError(f"失败单元：{cell['cell_id']}")
        if manifest["provenance"]["executable_sha256"] != expected_hash:
            raise RuntimeError(f"哈希不一致：{cell['cell_id']}")
        obs = manifest["observables"]
        rows.append(
            {
                "cell_id": cell["cell_id"],
                "temperature": float(manifest["params"]["temperature"]),
                "L": int(manifest["params"]["L"]),
                "binder": float(obs["binder_staggered"]["value"]),
                "binder_error": float(obs["binder_staggered"]["error"]),
                "energy": float(obs["energy_per_site"]["value"]),
                "energy_error": float(obs["energy_per_site"]["error"]),
                "lag1_energy": float(manifest["diagnostics"]["energy_bin_lag1"]),
                "lag1_m2": float(manifest["diagnostics"]["staggered_m2_bin_lag1"]),
            }
        )
    return sorted(rows, key=lambda row: (float(row["temperature"]), int(row["L"])))


def all_crossings(rows: list[dict[str, float | str]]) -> list[dict[str, float | int | list[float]]]:
    sizes = sorted({int(row["L"]) for row in rows})
    crossings = []
    for small, large in zip(sizes, sizes[1:]):
        a = {float(row["temperature"]): row for row in rows if int(row["L"]) == small}
        b = {float(row["temperature"]): row for row in rows if int(row["L"]) == large}
        temperatures = sorted(set(a) & set(b))
        differences = []
        for temperature in temperatures:
            delta = float(b[temperature]["binder"]) - float(a[temperature]["binder"])
            sigma = math.hypot(float(b[temperature]["binder_error"]), float(a[temperature]["binder_error"]))
            differences.append((temperature, delta, sigma))
        roots = []
        for left, right in zip(differences, differences[1:]):
            t1, d1, s1 = left
            t2, d2, s2 = right
            if d1 * d2 > 0.0:
                continue
            denominator = d2 - d1
            if denominator == 0.0:
                continue
            root = (t1 * d2 - t2 * d1) / denominator
            derivative1 = d2 * (t1 - t2) / denominator**2
            derivative2 = d1 * (t2 - t1) / denominator**2
            error = math.hypot(derivative1 * s1, derivative2 * s2)
            roots.append({"temperature": root, "error": error, "bracket": [t1, t2]})
        crossings.append({"L_small": small, "L_large": large, "roots": roots})
    return crossings


def write_csv(rows: list[dict[str, float | str]], output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_svg(rows: list[dict[str, float | str]], output: Path) -> None:
    width, height = 900, 600
    left, right, top, bottom = 96, 42, 42, 82
    plot_width, plot_height = width - left - right, height - top - bottom
    temperatures = [float(row["temperature"]) for row in rows]
    binders = [float(row["binder"]) for row in rows]
    xmin, xmax = min(temperatures) - 0.003, max(temperatures) + 0.003
    ymin, ymax = min(binders) - 0.06, min(0.68, max(binders) + 0.025)
    sizes = sorted({int(row["L"]) for row in rows})
    palette = ["#2166ac", "#4393c3", "#4d9221", "#d6604d", "#b2182b", "#762a83", "#e6ab02"]

    def x(value: float) -> float:
        return left + (value - xmin) * plot_width / (xmax - xmin)

    def y(value: float) -> float:
        return top + (ymax - value) * plot_height / (ymax - ymin)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.tick{font-size:14px}.label{font-size:18px}.legend{font-size:14px}</style>',
    ]
    for value in [0.4, 0.45, 0.5, 0.55, 0.6, 0.65]:
        if value < ymin or value > ymax:
            continue
        yy = y(value)
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{width-right}" y2="{yy:.2f}" stroke="#dddddd"/>')
        parts.append(f'<text class="tick" x="{left-12}" y="{yy+5:.2f}" text-anchor="end">{value:.2f}</text>')
    for value in sorted(set(temperatures)):
        xx = x(value)
        parts.append(f'<line x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{height-bottom}" stroke="#eeeeee"/>')
        parts.append(f'<text class="tick" x="{xx:.2f}" y="{height-bottom+28}" text-anchor="middle">{value:.3f}</text>')
    parts.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')
    exact = 2.0 / math.acosh(2.0 + math.sqrt(3.0))
    exact_x = x(exact)
    parts.append(f'<line x1="{exact_x:.2f}" y1="{top}" x2="{exact_x:.2f}" y2="{height-bottom}" stroke="#555" stroke-width="1.5" stroke-dasharray="7 5"/>')
    parts.append(f'<text class="legend" x="{exact_x+7:.2f}" y="{top+18}">exact Tc/J = {exact:.6f}</text>')

    for index, size in enumerate(sizes):
        color = palette[index % len(palette)]
        series = [row for row in rows if int(row["L"]) == size]
        points = " ".join(f'{x(float(row["temperature"])):.2f},{y(float(row["binder"])):.2f}' for row in series)
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>')
        for row in series:
            xx, yy = x(float(row["temperature"])), y(float(row["binder"]))
            y1 = y(float(row["binder"]) + float(row["binder_error"]))
            y2 = y(float(row["binder"]) - float(row["binder_error"]))
            parts.append(f'<line x1="{xx:.2f}" y1="{y1:.2f}" x2="{xx:.2f}" y2="{y2:.2f}" stroke="{color}"/>')
            parts.append(f'<line x1="{xx-3:.2f}" y1="{y1:.2f}" x2="{xx+3:.2f}" y2="{y1:.2f}" stroke="{color}"/>')
            parts.append(f'<line x1="{xx-3:.2f}" y1="{y2:.2f}" x2="{xx+3:.2f}" y2="{y2:.2f}" stroke="{color}"/>')
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.8" fill="white" stroke="{color}" stroke-width="2"/>')
        legend_x, legend_y = width - right - 92, top + 32 + 24 * index
        parts.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x+23}" y2="{legend_y}" stroke="{color}" stroke-width="2"/>')
        parts.append(f'<text class="legend" x="{legend_x+31}" y="{legend_y+5}">L={size}</text>')
    parts.append(f'<text class="label" x="{left+plot_width/2:.2f}" y="{height-24}" text-anchor="middle">Temperature T/J</text>')
    parts.append(f'<text class="label" transform="translate(28 {top+plot_height/2:.2f}) rotate(-90)" text-anchor="middle">Staggered Binder cumulant</text>')
    parts.append('</svg>')
    output.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = load_rows(args.run_dir)
    high_lag = [row for row in rows if abs(float(row["lag1_m2"])) >= 0.5]
    summary = {
        "cell_count": len(rows),
        "exact_temperature": 2.0 / math.acosh(2.0 + math.sqrt(3.0)),
        "adjacent_size_crossings": all_crossings(rows),
        "max_abs_staggered_m2_bin_lag1": max(abs(float(row["lag1_m2"])) for row in rows),
        "high_autocorrelation_cells": high_lag,
        "precision_status": "coarse_only" if high_lag else "statistically_usable",
    }
    write_csv(rows, args.run_dir / "scan_observables.csv")
    (args.run_dir / "scan_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_svg(rows, args.run_dir / "binder_curves.svg")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
