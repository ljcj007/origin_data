#!/usr/bin/env python3
"""计算周期格点图的最近邻 Bloch 邻接谱。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class GraphEdge:
    site1: int
    site2: int
    shift: tuple[int, int, int]
    bond_type: int


@dataclass(frozen=True)
class PeriodicGraph:
    basis_count: int
    edges: tuple[GraphEdge, ...]
    eta: tuple[int, ...]
    phase: tuple[int, int, int]


def read_graph(path: Path) -> PeriodicGraph:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    basis_count, edge_count, dimension = (int(value) for value in lines[0].split())
    if dimension != 3:
        raise ValueError("当前只支持三分量平移格式")
    if len(lines) != edge_count + 3:
        raise ValueError("格点图行数与文件头不一致")
    edges: list[GraphEdge] = []
    for line in lines[1 : edge_count + 1]:
        site1, site2, dx, dy, dz, bond_type = (int(value) for value in line.split())
        if not 1 <= site1 <= basis_count or not 1 <= site2 <= basis_count:
            raise ValueError("格点图含越界编号")
        edges.append(GraphEdge(site1, site2, (dx, dy, dz), bond_type))
    eta = tuple(int(value) for value in lines[edge_count + 1].split())
    phase = tuple(int(value) for value in lines[edge_count + 2].split())
    if len(eta) != basis_count or len(phase) != 3:
        raise ValueError("二分子格信息格式错误")
    return PeriodicGraph(basis_count, tuple(edges), eta, phase)  # type: ignore[arg-type]


def bloch_hamiltonian(graph: PeriodicGraph, kx: float, ky: float, hoppings: np.ndarray) -> np.ndarray:
    matrix = np.zeros((graph.basis_count, graph.basis_count), dtype=np.complex128)
    for edge in graph.edges:
        i = edge.site1 - 1
        j = edge.site2 - 1
        amplitude = hoppings[edge.bond_type - 1]
        phase = np.exp(1.0j * (kx * edge.shift[0] + ky * edge.shift[1]))
        matrix[i, j] += amplitude * phase
        matrix[j, i] += amplitude * phase.conjugate()
    return matrix


def reciprocal_path(
    path_kind: str, points_per_segment: int
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, str]]]:
    if path_kind == "square":
        vertices = [
            (np.array([0.0, 0.0]), "Γ"),
            (np.array([math.pi, 0.0]), "X"),
            (np.array([math.pi, math.pi]), "M"),
            (np.array([0.0, 0.0]), "Γ"),
        ]
    elif path_kind == "hexagonal":
        vertices = [
            (np.array([0.0, 0.0]), "Γ"),
            (np.array([math.pi, 0.0]), "M"),
            (np.array([2.0 * math.pi / 3.0, -2.0 * math.pi / 3.0]), "K"),
            (np.array([0.0, 0.0]), "Γ"),
        ]
    elif path_kind == "hexagonal120":
        vertices = [
            (np.array([0.0, 0.0]), "Γ"),
            (np.array([math.pi, 0.0]), "M"),
            (np.array([2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]), "K"),
            (np.array([0.0, 0.0]), "Γ"),
        ]
    else:
        raise ValueError(f"未知倒空间路径：{path_kind}")
    points: list[np.ndarray] = []
    distances: list[float] = []
    ticks: list[tuple[int, str]] = [(0, vertices[0][1])]
    distance = 0.0
    for segment in range(len(vertices) - 1):
        start, _ = vertices[segment]
        end, end_label = vertices[segment + 1]
        for step in range(points_per_segment):
            fraction = step / points_per_segment
            point = (1.0 - fraction) * start + fraction * end
            if points:
                distance += float(np.linalg.norm(point - points[-1]))
            points.append(point)
            distances.append(distance)
        distance += float(np.linalg.norm(end - points[-1]))
        points.append(end)
        distances.append(distance)
        ticks.append((len(points) - 1, end_label))
    return np.asarray(points), np.asarray(distances), ticks


def sample_grid(graph: PeriodicGraph, hoppings: np.ndarray, grid_size: int) -> np.ndarray:
    values: list[np.ndarray] = []
    grid = np.linspace(-math.pi, math.pi, grid_size, endpoint=False)
    for kx in grid:
        for ky in grid:
            values.append(np.linalg.eigvalsh(bloch_hamiltonian(graph, kx, ky, hoppings)))
    return np.asarray(values)


def detect_flat_energies(
    graph: PeriodicGraph,
    hoppings: np.ndarray,
    grid_values: np.ndarray,
    tolerance: float,
) -> list[dict[str, float | int]]:
    gamma_values = np.linalg.eigvalsh(bloch_hamiltonian(graph, 0.0, 0.0, hoppings))
    scale = max(1.0, float(np.max(np.abs(grid_values))))
    absolute_tolerance = tolerance * scale
    candidates: list[dict[str, float | int]] = []
    for energy in gamma_values:
        if any(abs(float(item["energy"]) - float(energy)) <= absolute_tolerance for item in candidates):
            continue
        nearest = np.min(np.abs(grid_values - energy), axis=1)
        residual = float(np.max(nearest))
        if residual <= absolute_tolerance:
            multiplicity = int(np.min(np.sum(np.abs(grid_values - energy) <= absolute_tolerance, axis=1)))
            candidates.append(
                {
                    "energy": float(energy),
                    "minimum_multiplicity": multiplicity,
                    "maximum_grid_residual": residual,
                }
            )
    return candidates


def write_path_csv(path: Path, distances: np.ndarray, points: np.ndarray, bands: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path_distance", "kx", "ky", *[f"band_{i + 1}" for i in range(bands.shape[1])]])
        for distance, point, energies in zip(distances, points, bands):
            writer.writerow([distance, point[0], point[1], *energies])


def svg_polyline(points: list[tuple[float, float]], color: str) -> str:
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="2"/>'


def write_svg(
    path: Path,
    distances: np.ndarray,
    bands: np.ndarray,
    ticks: list[tuple[int, str]],
    flat_energies: list[dict[str, float | int]],
) -> None:
    width, height = 900, 600
    left, right, top, bottom = 90, 30, 45, 75
    plot_width = width - left - right
    plot_height = height - top - bottom
    x_min, x_max = float(distances[0]), float(distances[-1])
    y_min, y_max = float(np.min(bands)), float(np.max(bands))
    padding = 0.06 * max(1.0, y_max - y_min)
    y_min -= padding
    y_max += padding

    def x_map(value: float) -> float:
        return left + plot_width * (value - x_min) / (x_max - x_min)

    def y_map(value: float) -> float:
        return top + plot_height * (y_max - value) / (y_max - y_min)

    colors = ["#2855a6", "#c23b33", "#2d8a55", "#8b4ba5", "#d17b18", "#168b92"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>',
        '<text x="450" y="28" text-anchor="middle" font-family="Arial" font-size="20">Nearest-neighbor Bloch spectrum</text>',
        f'<text x="22" y="{top + plot_height/2}" transform="rotate(-90 22 {top + plot_height/2})" text-anchor="middle" font-family="Arial" font-size="16">Eigenvalue / t</text>',
    ]
    for index in range(bands.shape[1]):
        line = [(x_map(float(x)), y_map(float(y))) for x, y in zip(distances, bands[:, index])]
        parts.append(svg_polyline(line, colors[index % len(colors)]))
    for tick_index, label in ticks:
        x = x_map(float(distances[tick_index]))
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#bbb"/>')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_height + 28}" text-anchor="middle" font-family="Arial" font-size="16">{label}</text>')
    for energy in flat_energies:
        y = y_map(float(energy["energy"]))
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#111" stroke-dasharray="6 5" opacity="0.55"/>')
    for fraction in np.linspace(0.0, 1.0, 7):
        value = y_min + fraction * (y_max - y_min)
        y = y_map(value)
        parts.append(f'<line x1="{left - 5}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#222"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-family="Arial" font-size="13">{value:.2f}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path)
    parser.add_argument("--hoppings", nargs="+", type=float, default=[1.0])
    parser.add_argument("--grid", type=int, default=41)
    parser.add_argument("--path-points", type=int, default=100)
    parser.add_argument("--path", choices=("square", "hexagonal", "hexagonal120"), default="square")
    parser.add_argument("--flat-tolerance", type=float, default=1.0e-9)
    parser.add_argument("--output-dir", type=Path, default=Path("results/bands"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    graph = read_graph(args.graph)
    maximum_type = max(edge.bond_type for edge in graph.edges)
    if len(args.hoppings) != maximum_type:
        raise SystemExit(f"需要为 {maximum_type} 类键分别给出 hopping")
    if args.grid < 3 or args.path_points < 2:
        raise SystemExit("网格和路径点数过小")

    hoppings = np.asarray(args.hoppings, dtype=float)
    points, distances, ticks = reciprocal_path(args.path, args.path_points)
    bands = np.asarray(
        [np.linalg.eigvalsh(bloch_hamiltonian(graph, point[0], point[1], hoppings)) for point in points]
    )
    grid_values = sample_grid(graph, hoppings, args.grid)
    flat_energies = detect_flat_energies(graph, hoppings, grid_values, args.flat_tolerance)
    sorted_widths = np.max(grid_values, axis=0) - np.min(grid_values, axis=0)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "band_path.csv"
    svg_path = args.output_dir / "band_structure.svg"
    report_path = args.output_dir / "flat_band_report.json"
    write_path_csv(csv_path, distances, points, bands)
    write_svg(svg_path, distances, bands, ticks, flat_energies)
    report = {
        "format": "periodic-graph-band-report-v1",
        "graph": str(args.graph.resolve()),
        "basis_count": graph.basis_count,
        "hoppings": args.hoppings,
        "grid": args.grid,
        "flat_tolerance": args.flat_tolerance,
        "sorted_band_widths": sorted_widths.tolist(),
        "flat_energies": flat_energies,
        "path": args.path,
        "path_note": (
            "路径使用晶胞倒空间相位坐标；hexagonal120 对应面内基矢夹角 120 度。"
            if args.path == "hexagonal120"
            else "路径使用晶胞倒空间相位坐标。"
        ),
    }
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"基底格点数：{graph.basis_count}")
    print(f"排序能带宽度：{' '.join(f'{value:.10g}' for value in sorted_widths)}")
    if flat_energies:
        for item in flat_energies:
            print(
                f"平带：E={item['energy']:.12g}, 最小简并度={item['minimum_multiplicity']}, "
                f"最大残差={item['maximum_grid_residual']:.3e}"
            )
    else:
        print("在给定容差和网格上未检测到严格平带")
    print(f"报告：{report_path}")
    print(f"图：{svg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
