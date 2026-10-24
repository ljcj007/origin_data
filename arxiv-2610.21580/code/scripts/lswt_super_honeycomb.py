#!/usr/bin/env python3
"""计算边心蜂窝反铁磁 Heisenberg 模型的线性自旋波谱。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Edge:
    site1: int
    site2: int
    shift: tuple[int, int, int]


@dataclass(frozen=True)
class Graph:
    basis_count: int
    edges: tuple[Edge, ...]
    eta: tuple[int, ...]


def read_graph(path: Path) -> Graph:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    basis_count, edge_count, dimension = (int(value) for value in lines[0].split())
    if dimension != 3 or len(lines) != edge_count + 3:
        raise ValueError("格点图格式错误")
    edges: list[Edge] = []
    for line in lines[1 : edge_count + 1]:
        site1, site2, dx, dy, dz, _ = (int(value) for value in line.split())
        edges.append(Edge(site1 - 1, site2 - 1, (dx, dy, dz)))
    eta = tuple(int(value) for value in lines[edge_count + 1].split())
    if len(eta) != basis_count or any(value not in (-1, 1) for value in eta):
        raise ValueError("需要已标记的二分图")
    return Graph(basis_count, tuple(edges), eta)


def orient_bonds(graph: Graph) -> tuple[list[int], list[int], list[tuple[int, int, tuple[int, int]]]]:
    a_sites = [index for index, value in enumerate(graph.eta) if value == 1]
    b_sites = [index for index, value in enumerate(graph.eta) if value == -1]
    a_index = {site: index for index, site in enumerate(a_sites)}
    b_index = {site: index for index, site in enumerate(b_sites)}
    bonds: list[tuple[int, int, tuple[int, int]]] = []
    for edge in graph.edges:
        if graph.eta[edge.site1] * graph.eta[edge.site2] != -1:
            raise ValueError("发现同子格键，当前 LSWT 推导不适用")
        dx, dy = edge.shift[:2]
        if graph.eta[edge.site1] == 1:
            bonds.append((a_index[edge.site1], b_index[edge.site2], (dx, dy)))
        else:
            bonds.append((a_index[edge.site2], b_index[edge.site1], (-dx, -dy)))
    return a_sites, b_sites, bonds


def coordination(
    n_a: int, n_b: int, bonds: list[tuple[int, int, tuple[int, int]]]
) -> tuple[np.ndarray, np.ndarray]:
    z_a = np.zeros(n_a, dtype=float)
    z_b = np.zeros(n_b, dtype=float)
    for a_index, b_index, _ in bonds:
        z_a[a_index] += 1.0
        z_b[b_index] += 1.0
    return z_a, z_b


def coupling_matrix(
    n_a: int,
    n_b: int,
    bonds: list[tuple[int, int, tuple[int, int]]],
    kx: float,
    ky: float,
) -> np.ndarray:
    matrix = np.zeros((n_a, n_b), dtype=np.complex128)
    for a_index, b_index, (dx, dy) in bonds:
        matrix[a_index, b_index] += np.exp(1.0j * (kx * dx + ky * dy))
    return matrix


def branch_energies(
    coupling: np.ndarray, z_a: float, z_b: float, spin: float, exchange: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sigma2 = np.linalg.svd(coupling, compute_uv=False) ** 2
    sigma2 = np.sort(sigma2)[::-1]
    delta = z_a - z_b
    radicands = (z_a + z_b) ** 2 - 4.0 * sigma2
    if float(np.min(radicands)) < -1.0e-10:
        raise ValueError("出现不稳定的虚频")
    roots = np.sqrt(np.maximum(radicands, 0.0))
    minus = 2.0 * np.maximum(z_a * z_b - sigma2, 0.0) / (roots + delta)
    plus = 0.5 * (roots + delta)
    scale = spin * exchange
    flat = np.full(coupling.shape[1] - coupling.shape[0], z_b)
    labeled = np.concatenate((minus, plus, flat)) * scale
    return labeled, np.sort(labeled), sigma2


def dynamical_eigenvalues(coupling: np.ndarray, z_a: float, z_b: float) -> np.ndarray:
    n_a, n_b = coupling.shape
    matrix = np.block(
        [
            [z_a * np.eye(n_a), -coupling],
            [coupling.conjugate().T, -z_b * np.eye(n_b)],
        ]
    )
    values = np.linalg.eigvals(matrix)
    if float(np.max(np.abs(values.imag))) > 1.0e-9:
        raise ValueError("动力学矩阵含显著虚频")
    return np.sort(values.real)


def reciprocal_path(points_per_segment: int) -> tuple[np.ndarray, np.ndarray, list[tuple[int, str]]]:
    vertices = [
        (np.array([0.0, 0.0]), "Γ"),
        (np.array([math.pi, 0.0]), "M"),
        (np.array([2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]), "K"),
        (np.array([0.0, 0.0]), "Γ"),
    ]
    points: list[np.ndarray] = []
    distances: list[float] = []
    ticks: list[tuple[int, str]] = [(0, "Γ")]
    distance = 0.0
    for segment in range(3):
        start, _ = vertices[segment]
        end, label = vertices[segment + 1]
        for step in range(points_per_segment):
            point = start + (end - start) * step / points_per_segment
            if points:
                distance += float(np.linalg.norm(point - points[-1]))
            points.append(point)
            distances.append(distance)
        distance += float(np.linalg.norm(end - points[-1]))
        points.append(end)
        distances.append(distance)
        ticks.append((len(points) - 1, label))
    return np.asarray(points), np.asarray(distances), ticks


def write_path_csv(
    path: Path,
    points: np.ndarray,
    distances: np.ndarray,
    labeled: np.ndarray,
    sorted_bands: np.ndarray,
    sigma2: np.ndarray,
) -> None:
    header = [
        "path_distance",
        "kx",
        "ky",
        "omega_sigma_max_minus",
        "omega_sigma_min_minus",
        "omega_sigma_max_plus",
        "omega_sigma_min_plus",
        "omega_flat",
        "omega_sorted_1",
        "omega_sorted_2",
        "omega_sorted_3",
        "omega_sorted_4",
        "omega_sorted_5",
        "sigma_max_squared",
        "sigma_min_squared",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for distance, point, row_labeled, row_sorted, row_sigma2 in zip(
            distances, points, labeled, sorted_bands, sigma2
        ):
            writer.writerow([distance, *point, *row_labeled, *row_sorted, *row_sigma2])


def svg_polyline(points: list[tuple[float, float]], color: str, dash: str = "") -> str:
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<polyline points="{coordinates}" fill="none" stroke="{color}" '
        f'stroke-width="2.3"{dash_attr}/>'
    )


def write_svg(
    path: Path,
    distances: np.ndarray,
    labeled: np.ndarray,
    ticks: list[tuple[int, str]],
    spin: float,
    exchange: float,
) -> None:
    width, height = 900, 610
    left, right, top, bottom = 92, 35, 48, 78
    plot_width = width - left - right
    plot_height = height - top - bottom
    x_min, x_max = float(distances[0]), float(distances[-1])
    y_min, y_max = 0.0, float(np.max(labeled)) * 1.06

    def x_map(value: float) -> float:
        return left + plot_width * (value - x_min) / (x_max - x_min)

    def y_map(value: float) -> float:
        return top + plot_height * (y_max - value) / (y_max - y_min)

    colors = ["#005f73", "#ca6702", "#0a9396", "#bb3e03", "#5a189a"]
    labels = ["σmax −", "σmin −", "σmax +", "σmin +", "flat"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>',
        '<text x="450" y="29" text-anchor="middle" font-family="Arial" font-size="20">Super-honeycomb LSWT magnon bands</text>',
        f'<text x="24" y="{top + plot_height / 2}" transform="rotate(-90 24 {top + plot_height / 2})" text-anchor="middle" font-family="Arial" font-size="16">ω / J (S={spin:g})</text>',
    ]
    for band in range(labeled.shape[1]):
        line = [(x_map(float(x)), y_map(float(y))) for x, y in zip(distances, labeled[:, band])]
        parts.append(svg_polyline(line, colors[band], "7 5" if band == 4 else ""))
    for tick_index, label in ticks:
        x_value = x_map(float(distances[tick_index]))
        parts.append(
            f'<line x1="{x_value:.2f}" y1="{top}" x2="{x_value:.2f}" y2="{top + plot_height}" stroke="#c8c8c8"/>'
        )
        parts.append(
            f'<text x="{x_value:.2f}" y="{top + plot_height + 29}" text-anchor="middle" font-family="Arial" font-size="17">{label}</text>'
        )
    for fraction in np.linspace(0.0, 1.0, 7):
        value = y_min + fraction * (y_max - y_min)
        y_value = y_map(value)
        parts.append(f'<line x1="{left - 5}" y1="{y_value:.2f}" x2="{left}" y2="{y_value:.2f}" stroke="#222"/>')
        parts.append(
            f'<text x="{left - 10}" y="{y_value + 5:.2f}" text-anchor="end" font-family="Arial" font-size="13">{value / exchange:.2f}</text>'
        )
    legend_x, legend_y = left + 18, top + 22
    for index, (color, label) in enumerate(zip(colors, labels)):
        y_value = legend_y + 22 * index
        dash = ' stroke-dasharray="7 5"' if index == 4 else ""
        parts.append(f'<line x1="{legend_x}" y1="{y_value}" x2="{legend_x + 31}" y2="{y_value}" stroke="{color}" stroke-width="2.3"{dash}/>')
        parts.append(f'<text x="{legend_x + 39}" y="{y_value + 5}" font-family="Arial" font-size="13">{label}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_small_k_csv(
    path: Path,
    n_a: int,
    n_b: int,
    bonds: list[tuple[int, int, tuple[int, int]]],
    z_a: float,
    z_b: float,
    spin: float,
    exchange: float,
) -> tuple[float, float]:
    q_values = np.geomspace(1.0e-4, 8.0e-2, 48)
    momenta: list[float] = []
    energies: list[float] = []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["kx", "ky", "k_cartesian_a", "omega_acoustic_over_J"])
        for q_value in q_values:
            coupling = coupling_matrix(n_a, n_b, bonds, q_value, 0.0)
            labeled, _, _ = branch_energies(coupling, z_a, z_b, spin, exchange)
            momentum = 2.0 * q_value / math.sqrt(3.0)
            energy = float(labeled[0] / exchange)
            writer.writerow([q_value, 0.0, momentum, energy])
            momenta.append(momentum)
            energies.append(energy)
    log_fit = np.polyfit(np.log(momenta), np.log(energies), 1)
    coefficient = float(np.median(np.asarray(energies[:16]) / np.square(momenta[:16])))
    return float(log_fit[0]), coefficient


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=Path("graphs/super_honeycomb_poscar.dat"))
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument("--exchange", type=float, default=1.0)
    parser.add_argument("--grid", type=int, default=121)
    parser.add_argument("--path-points", type=int, default=160)
    parser.add_argument("--output-dir", type=Path, default=Path("results/lswt_super_honeycomb_v1"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.spin <= 0.0 or args.exchange <= 0.0 or args.grid < 5 or args.path_points < 4:
        raise SystemExit("参数超出允许范围")
    graph = read_graph(args.graph)
    a_sites, b_sites, bonds = orient_bonds(graph)
    if len(a_sites) >= len(b_sites):
        raise SystemExit("约定要求 eta=+1 为较小子格")
    z_a_values, z_b_values = coordination(len(a_sites), len(b_sites), bonds)
    if np.ptp(z_a_values) > 0.0 or np.ptp(z_b_values) > 0.0:
        raise SystemExit("当前解析分支要求各子格内部配位数相同")
    z_a, z_b = float(z_a_values[0]), float(z_b_values[0])

    points, distances, ticks = reciprocal_path(args.path_points)
    labeled_rows: list[np.ndarray] = []
    sorted_rows: list[np.ndarray] = []
    sigma2_rows: list[np.ndarray] = []
    dynamic_residual = 0.0
    for point in points:
        coupling = coupling_matrix(len(a_sites), len(b_sites), bonds, point[0], point[1])
        labeled, sorted_bands, sigma2 = branch_energies(coupling, z_a, z_b, args.spin, args.exchange)
        signed_expected = np.concatenate(
            (
                0.5 * (z_a - z_b - np.sqrt((z_a + z_b) ** 2 - 4.0 * sigma2)),
                0.5 * (z_a - z_b + np.sqrt((z_a + z_b) ** 2 - 4.0 * sigma2)),
                np.full(len(b_sites) - len(a_sites), -z_b),
            )
        )
        direct = dynamical_eigenvalues(coupling, z_a, z_b)
        dynamic_residual = max(dynamic_residual, float(np.max(np.abs(direct - np.sort(signed_expected)))))
        labeled_rows.append(labeled)
        sorted_rows.append(sorted_bands)
        sigma2_rows.append(sigma2)
    labeled_array = np.asarray(labeled_rows)
    sorted_array = np.asarray(sorted_rows)
    sigma2_array = np.asarray(sigma2_rows)

    flat_residual = 0.0
    minimum_frequency = math.inf
    grid_values = 2.0 * math.pi * (np.arange(args.grid) - args.grid // 2) / args.grid
    for kx in grid_values:
        for ky in grid_values:
            coupling = coupling_matrix(len(a_sites), len(b_sites), bonds, kx, ky)
            labeled, _, _ = branch_energies(coupling, z_a, z_b, args.spin, args.exchange)
            direct = dynamical_eigenvalues(coupling, z_a, z_b)
            flat_residual = max(
                flat_residual,
                float(np.min(np.abs(direct + z_b))) * args.spin * args.exchange,
            )
            minimum_frequency = min(minimum_frequency, float(np.min(labeled)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_path_csv(
        args.output_dir / "magnon_path.csv",
        points,
        distances,
        labeled_array,
        sorted_array,
        sigma2_array,
    )
    write_svg(args.output_dir / "magnon_bands.svg", distances, labeled_array, ticks, args.spin, args.exchange)
    exponent, coefficient = write_small_k_csv(
        args.output_dir / "small_k_dispersion.csv",
        len(a_sites),
        len(b_sites),
        bonds,
        z_a,
        z_b,
        args.spin,
        args.exchange,
    )

    special_points: dict[str, list[float]] = {}
    for label, point in {
        "Gamma": (0.0, 0.0),
        "M": (math.pi, 0.0),
        "K": (2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0),
    }.items():
        coupling = coupling_matrix(len(a_sites), len(b_sites), bonds, *point)
        _, sorted_bands, _ = branch_energies(coupling, z_a, z_b, args.spin, args.exchange)
        special_points[label] = [float(value / args.exchange) for value in sorted_bands]

    expected_coefficient = args.spin / 4.0
    summary = {
        "format": "super-honeycomb-lswt-v1",
        "model": "H = J sum_<ij> S_i dot S_j",
        "convention": {"J": args.exchange, "S": args.spin, "J_positive_is_antiferromagnetic": True},
        "input": {"graph": str(args.graph.resolve()), "graph_sha256": sha256(args.graph)},
        "lattice": {
            "basis_count": graph.basis_count,
            "A_sites_one_based": [value + 1 for value in a_sites],
            "B_sites_one_based": [value + 1 for value in b_sites],
            "A_coordination": z_a,
            "B_coordination": z_b,
            "classical_total_spin_per_cell": (len(b_sites) - len(a_sites)) * args.spin,
        },
        "bands": {
            "count": graph.basis_count,
            "flat_band_count": len(b_sites) - len(a_sites),
            "flat_energy_over_J": z_b * args.spin,
            "flat_grid_residual_over_J": flat_residual / args.exchange,
            "special_point_frequencies_over_J": special_points,
        },
        "goldstone": {
            "count": 1,
            "type": "type-B quadratic",
            "linear_mode_present": False,
            "fitted_dynamic_exponent": exponent,
            "quadratic_coefficient_over_J_in_unit_primitive_cartesian_coordinates": coefficient,
            "analytic_quadratic_coefficient": expected_coefficient,
        },
        "verification": {
            "max_reduced_dynamic_matrix_residual": dynamic_residual,
            "minimum_grid_frequency_over_J": minimum_frequency / args.exchange,
            "quadratic_coefficient_relative_error": abs(coefficient - expected_coefficient) / expected_coefficient,
        },
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"五条磁振子带，平带 ω/J={z_b * args.spin:.12g}")
    print(f"Goldstone 指数 z={exponent:.8f}，二次项系数={coefficient:.8f}")
    print(f"动力学矩阵最大残差={dynamic_residual:.3e}")
    print(f"结果目录：{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
