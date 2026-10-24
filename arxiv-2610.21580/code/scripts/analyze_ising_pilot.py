#!/usr/bin/env python3
"""汇总 Ising 预扫描并绘制 Binder 交叉图。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def load_rows(run_dir: Path) -> list[dict[str, float]]:
    spec = json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))
    expected_hash = spec["provenance"]["executable_sha256"]
    rows = []
    for cell in spec["cells"]:
        path = run_dir / "cells" / cell["cell_id"] / "manifest.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        if not data["status"]["success"]:
            raise RuntimeError(f"单元失败：{cell['cell_id']}")
        if data["provenance"]["executable_sha256"] != expected_hash:
            raise RuntimeError(f"可执行文件哈希不一致：{cell['cell_id']}")
        obs = data["observables"]
        rows.append(
            {
                "temperature": float(data["params"]["temperature"]),
                "L": float(data["params"]["L"]),
                "binder": float(obs["binder_staggered"]["value"]),
                "binder_error": float(obs["binder_staggered"]["error"]),
                "energy": float(obs["energy_per_site"]["value"]),
                "energy_error": float(obs["energy_per_site"]["error"]),
                "acceptance": float(obs["acceptance_rate"]["value"]),
            }
        )
    return sorted(rows, key=lambda row: (row["temperature"], row["L"]))


def estimate_crossing(rows: list[dict[str, float]]) -> tuple[float, float, tuple[float, float]]:
    sizes = sorted({int(row["L"]) for row in rows})
    if len(sizes) != 2:
        raise RuntimeError("预扫描交叉估计要求恰好两个尺寸")
    by_size = {
        size: {row["temperature"]: row for row in rows if int(row["L"]) == size}
        for size in sizes
    }
    temperatures = sorted(set(by_size[sizes[0]]) & set(by_size[sizes[1]]))
    differences = []
    for temperature in temperatures:
        small = by_size[sizes[0]][temperature]
        large = by_size[sizes[1]][temperature]
        delta = large["binder"] - small["binder"]
        sigma = math.hypot(large["binder_error"], small["binder_error"])
        differences.append((temperature, delta, sigma))

    for left, right in zip(differences, differences[1:]):
        t1, d1, s1 = left
        t2, d2, s2 = right
        if d1 == 0.0 or d1 * d2 <= 0.0:
            denominator = d2 - d1
            crossing = (t1 * d2 - t2 * d1) / denominator
            derivative1 = d2 * (t1 - t2) / denominator**2
            derivative2 = d1 * (t2 - t1) / denominator**2
            error = math.hypot(derivative1 * s1, derivative2 * s2)
            return crossing, error, (t1, t2)
    raise RuntimeError("扫描区间内没有 Binder 交叉")


def write_csv(rows: list[dict[str, float]], output: Path) -> None:
    fields = list(rows[0])
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_svg(rows: list[dict[str, float]], crossing: float, crossing_error: float, output: Path) -> None:
    width, height = 820, 560
    left, right, top, bottom = 92, 36, 42, 78
    plot_width = width - left - right
    plot_height = height - top - bottom
    xmin, xmax = 0.89, 1.11
    ymin, ymax = 0.18, 0.68

    def x(value: float) -> float:
        return left + (value - xmin) * plot_width / (xmax - xmin)

    def y(value: float) -> float:
        return top + (ymax - value) * plot_height / (ymax - ymin)

    exact = 2.0 / math.acosh(2.0 + math.sqrt(3.0))
    colors = {8: "#2166ac", 12: "#b2182b"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.tick{font-size:14px}.label{font-size:18px}.legend{font-size:15px}</style>',
    ]
    for value in [0.2, 0.3, 0.4, 0.5, 0.6]:
        yy = y(value)
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{width-right}" y2="{yy:.2f}" stroke="#dddddd"/>')
        parts.append(f'<text class="tick" x="{left-12}" y="{yy+5:.2f}" text-anchor="end">{value:.1f}</text>')
    for value in [0.90, 0.95, 1.00, 1.05, 1.10]:
        xx = x(value)
        parts.append(f'<line x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{height-bottom}" stroke="#eeeeee"/>')
        parts.append(f'<text class="tick" x="{xx:.2f}" y="{height-bottom+28}" text-anchor="middle">{value:.2f}</text>')
    parts.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')

    exact_x = x(exact)
    parts.append(f'<line x1="{exact_x:.2f}" y1="{top}" x2="{exact_x:.2f}" y2="{height-bottom}" stroke="#555" stroke-width="1.5" stroke-dasharray="7 5"/>')
    parts.append(f'<text class="legend" x="{exact_x+7:.2f}" y="{top+20}">exact Tc/J = {exact:.5f}</text>')

    for size in sorted(colors):
        series = [row for row in rows if int(row["L"]) == size]
        points = " ".join(f'{x(row["temperature"]):.2f},{y(row["binder"]):.2f}' for row in series)
        color = colors[size]
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        for row in series:
            xx, yy = x(row["temperature"]), y(row["binder"])
            y1, y2 = y(row["binder"] + row["binder_error"]), y(row["binder"] - row["binder_error"])
            parts.append(f'<line x1="{xx:.2f}" y1="{y1:.2f}" x2="{xx:.2f}" y2="{y2:.2f}" stroke="{color}"/>')
            parts.append(f'<line x1="{xx-4:.2f}" y1="{y1:.2f}" x2="{xx+4:.2f}" y2="{y1:.2f}" stroke="{color}"/>')
            parts.append(f'<line x1="{xx-4:.2f}" y1="{y2:.2f}" x2="{xx+4:.2f}" y2="{y2:.2f}" stroke="{color}"/>')
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="4.5" fill="white" stroke="{color}" stroke-width="2"/>')

    legend_x, legend_y = width - right - 165, top + 48
    for offset, size in enumerate(sorted(colors)):
        yy = legend_y + 26 * offset
        parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x+30}" y2="{yy}" stroke="{colors[size]}" stroke-width="2.2"/>')
        parts.append(f'<circle cx="{legend_x+15}" cy="{yy}" r="4" fill="white" stroke="{colors[size]}" stroke-width="2"/>')
        parts.append(f'<text class="legend" x="{legend_x+40}" y="{yy+5}">L = {size}</text>')

    parts.append(f'<text class="label" x="{left+plot_width/2:.2f}" y="{height-22}" text-anchor="middle">Temperature T/J</text>')
    parts.append(f'<text class="label" transform="translate(27 {top+plot_height/2:.2f}) rotate(-90)" text-anchor="middle">Staggered Binder cumulant</text>')
    parts.append(f'<text class="legend" x="{left+10}" y="{height-bottom-16}">linear crossing: {crossing:.5f} ± {crossing_error:.5f}</text>')
    parts.append('</svg>')
    output.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = load_rows(args.run_dir)
    crossing, error, bracket = estimate_crossing(rows)
    exact = 2.0 / math.acosh(2.0 + math.sqrt(3.0))
    write_csv(rows, args.run_dir / "pilot_observables.csv")
    summary = {
        "binder_crossing_linear": crossing,
        "binder_crossing_error_propagated": error,
        "crossing_bracket": list(bracket),
        "exact_temperature": exact,
        "deviation_from_exact": crossing - exact,
        "note": "仅为两尺寸预扫描；误差未包含有限尺寸漂移。",
    }
    (args.run_dir / "pilot_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_svg(rows, crossing, error, args.run_dir / "binder_crossing.svg")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
