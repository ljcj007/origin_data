#!/usr/bin/env python3
"""计数边心蜂窝有限环面的零模、CLS 与非收缩线态。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from band_graph import PeriodicGraph, bloch_hamiltonian, read_graph


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH = ROOT / "graphs" / "super_honeycomb_poscar.dat"
DEFAULT_OUTPUT = ROOT / "results" / "flatband_cls_topology_v1"
RANK_PRIME = 1_000_003


@dataclass(frozen=True)
class FiniteGraph:
    size: int
    vertex_sites: tuple[int, ...]
    center_sites: tuple[int, ...]
    incidence: np.ndarray

    @property
    def cells(self) -> int:
        return self.size * self.size


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def cell_index(x: int, y: int, size: int) -> int:
    return (x % size) * size + (y % size)


def site_degrees(graph: PeriodicGraph) -> np.ndarray:
    degrees = np.zeros(graph.basis_count, dtype=int)
    for edge in graph.edges:
        degrees[edge.site1 - 1] += 1
        degrees[edge.site2 - 1] += 1
    return degrees


def build_finite_graph(graph: PeriodicGraph, size: int) -> FiniteGraph:
    """构造顶点—边心关联矩阵。"""
    if size < 2:
        raise ValueError("L 必须至少为 2")
    degrees = site_degrees(graph)
    vertex_sites = tuple(int(i + 1) for i in np.flatnonzero(degrees == 3))
    center_sites = tuple(int(i + 1) for i in np.flatnonzero(degrees == 2))
    if len(vertex_sites) != 2 or len(center_sites) != 3:
        raise ValueError("目标图必须含两个三配位顶点和三个二配位边心点")
    if any(graph.eta[site - 1] != 1 for site in vertex_sites) or any(
        graph.eta[site - 1] != -1 for site in center_sites
    ):
        raise ValueError("目标图的二分子格标签不符合约定")

    cells = size * size
    vertex_position = {site: index for index, site in enumerate(vertex_sites)}
    center_position = {site: index for index, site in enumerate(center_sites)}
    incidence = np.zeros((2 * cells, 3 * cells), dtype=np.int64)
    for x in range(size):
        for y in range(size):
            for edge in graph.edges:
                site1, site2 = edge.site1, edge.site2
                dx, dy, dz = edge.shift
                if dz != 0:
                    raise ValueError("当前计数只支持二维键")
                if site1 in vertex_position and site2 in center_position:
                    vertex_site, center_site = site1, site2
                    vertex_cell = cell_index(x, y, size)
                    center_cell = cell_index(x + dx, y + dy, size)
                elif site2 in vertex_position and site1 in center_position:
                    vertex_site, center_site = site2, site1
                    vertex_cell = cell_index(x + dx, y + dy, size)
                    center_cell = cell_index(x, y, size)
                else:
                    raise ValueError("发现非顶点—边心键")
                row = vertex_position[vertex_site] * cells + vertex_cell
                column = center_position[center_site] * cells + center_cell
                incidence[row, column] += 1
    if np.any(np.sum(incidence, axis=1) != 3) or np.any(np.sum(incidence, axis=0) != 2):
        raise ValueError("有限周期图的配位数不正确")
    return FiniteGraph(size, vertex_sites, center_sites, incidence)


def modular_rank(matrix: np.ndarray, prime: int = RANK_PRIME) -> int:
    """用有限域消元给出整数矩阵秩的下界证书。"""
    work = np.mod(matrix, prime).astype(np.int64, copy=True)
    rows, columns = work.shape
    pivot_row = 0
    for column in range(columns):
        candidates = np.flatnonzero(work[pivot_row:, column])
        if candidates.size == 0:
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            work[[pivot_row, selected]] = work[[selected, pivot_row]]
        inverse = pow(int(work[pivot_row, column]), -1, prime)
        work[pivot_row, column:] = (work[pivot_row, column:] * inverse) % prime
        if pivot_row + 1 < rows:
            factors = work[pivot_row + 1 :, column].copy()
            active = np.flatnonzero(factors)
            if active.size:
                targets = pivot_row + 1 + active
                work[targets, column:] = (
                    work[targets, column:]
                    - factors[active, None] * work[pivot_row, column:][None, :]
                ) % prime
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def center_index(finite: FiniteGraph, site: int, x: int, y: int) -> int:
    position = finite.center_sites.index(site)
    return position * finite.cells + cell_index(x, y, finite.size)


def vertex_index(finite: FiniteGraph, site: int, x: int, y: int) -> int:
    position = finite.vertex_sites.index(site)
    return position * finite.cells + cell_index(x, y, finite.size)


def construct_face_cls(finite: FiniteGraph) -> np.ndarray:
    """构造每个六边形上的交替振幅 CLS。"""
    if finite.center_sites != (3, 4, 5):
        raise ValueError("六边形公式要求目标图的 3、4、5 号边心编号")
    vectors = np.zeros((3 * finite.cells, finite.cells), dtype=np.int64)
    column = 0
    for x in range(finite.size):
        for y in range(finite.size):
            entries = (
                (3, x, y, 1),
                (4, x, y, -1),
                (5, x - 1, y, 1),
                (3, x - 1, y, -1),
                (4, x - 1, y - 1, 1),
                (5, x - 1, y - 1, -1),
            )
            for site, cell_x, cell_y, value in entries:
                vectors[center_index(finite, site, cell_x, cell_y), column] += value
            column += 1
    return vectors


def construct_line_states(finite: FiniteGraph) -> np.ndarray:
    """构造沿两条基本周期方向的非收缩线态。"""
    vectors = np.zeros((3 * finite.cells, 2), dtype=np.int64)
    for y in range(finite.size):
        vectors[center_index(finite, 3, 0, y), 0] = 1
        vectors[center_index(finite, 4, 0, y), 0] = -1
    for x in range(finite.size):
        vectors[center_index(finite, 4, x, 0), 1] = 1
        vectors[center_index(finite, 5, x, 0), 1] = -1
    return vectors


def construct_global_vertex_state(finite: FiniteGraph) -> np.ndarray:
    """构造只占据两个蜂窝顶点子格的全局零模。"""
    vector = np.zeros((2 * finite.cells, 1), dtype=np.int64)
    for x in range(finite.size):
        for y in range(finite.size):
            vector[vertex_index(finite, finite.vertex_sites[0], x, y), 0] = 1
            vector[vertex_index(finite, finite.vertex_sites[1], x, y), 0] = -1
    return vector


def full_adjacency(finite: FiniteGraph) -> np.ndarray:
    zeros_vertex = np.zeros((2 * finite.cells, 2 * finite.cells), dtype=np.int64)
    zeros_center = np.zeros((3 * finite.cells, 3 * finite.cells), dtype=np.int64)
    return np.block(
        [[zeros_vertex, finite.incidence], [finite.incidence.T, zeros_center]]
    )


def bloch_zero_count(graph: PeriodicGraph, size: int, tolerance: float) -> tuple[int, int, int, float]:
    total = 0
    gamma = 0
    nongamma: list[int] = []
    maximum_residual = 0.0
    for nx in range(size):
        for ny in range(size):
            kx = 2.0 * math.pi * nx / size
            ky = 2.0 * math.pi * ny / size
            values = np.linalg.eigvalsh(bloch_hamiltonian(graph, kx, ky, np.ones(1)))
            zero_values = np.abs(values[np.abs(values) <= tolerance])
            count = int(zero_values.size)
            total += count
            if zero_values.size:
                maximum_residual = max(maximum_residual, float(np.max(zero_values)))
            if nx == 0 and ny == 0:
                gamma = count
            else:
                nongamma.append(count)
    return total, gamma, min(nongamma), maximum_residual


def analyze_size(graph: PeriodicGraph, size: int, tolerance: float = 1.0e-10) -> dict[str, Any]:
    finite = build_finite_graph(graph, size)
    cells = finite.cells
    incidence_rank = modular_rank(finite.incidence)
    face_cls = construct_face_cls(finite)
    line_states = construct_line_states(finite)
    global_state = construct_global_vertex_state(finite)
    cls_rank = modular_rank(face_cls)
    edge_basis = np.column_stack((face_cls, line_states))
    edge_basis_rank = modular_rank(edge_basis)
    complete_basis = np.block(
        [
            [np.zeros((2 * cells, edge_basis.shape[1]), dtype=np.int64), global_state],
            [edge_basis, np.zeros((3 * cells, 1), dtype=np.int64)],
        ]
    )
    complete_rank = modular_rank(complete_basis)
    adjacency = full_adjacency(finite)
    direct_values = np.linalg.eigvalsh(adjacency.astype(float))
    direct_zero_count = int(np.count_nonzero(np.abs(direct_values) <= tolerance))
    direct_zero_residual = float(np.max(np.abs(direct_values[np.abs(direct_values) <= tolerance])))
    bloch_total, gamma_zero, nongamma_zero, bloch_residual = bloch_zero_count(
        graph, size, tolerance
    )

    cls_residual = int(np.max(np.abs(finite.incidence @ face_cls)))
    line_residual = int(np.max(np.abs(finite.incidence @ line_states)))
    global_residual = int(np.max(np.abs(finite.incidence.T @ global_state)))
    expected_incidence_rank = 2 * cells - 1
    expected_edge_nullity = cells + 1
    expected_full_nullity = cells + 2
    checks = {
        "incidence_rank": incidence_rank == expected_incidence_rank,
        "cls_rank": cls_rank == cells - 1,
        "edge_basis_complete": edge_basis_rank == expected_edge_nullity,
        "complete_basis": complete_rank == expected_full_nullity,
        "direct_spectrum": direct_zero_count == expected_full_nullity,
        "bloch_sum": bloch_total == expected_full_nullity,
        "gamma_multiplicity": gamma_zero == 3 and nongamma_zero == 1,
        "annihilation": cls_residual == line_residual == global_residual == 0,
        "single_cls_dependency": bool(np.all(np.sum(face_cls, axis=1) == 0)),
    }
    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise RuntimeError(f"L={size} 计数验证失败：{failed}")

    return {
        "L": size,
        "cells": cells,
        "sites": 5 * cells,
        "incidence_shape": list(finite.incidence.shape),
        "incidence_rank": incidence_rank,
        "edge_center_zero_dimension": 3 * cells - incidence_rank,
        "vertex_zero_dimension": 2 * cells - incidence_rank,
        "full_zero_dimension": 5 * cells - 2 * incidence_rank,
        "face_cls_count": cells,
        "face_cls_span_rank": cls_rank,
        "face_cls_dependency_count": cells - cls_rank,
        "noncontractible_line_count": 2,
        "cls_plus_lines_rank": edge_basis_rank,
        "global_vertex_state_count": 1,
        "complete_constructed_rank": complete_rank,
        "nominal_flat_band_states": cells,
        "gamma_touching_extra_zero_modes": expected_full_nullity - cells,
        "gamma_zero_multiplicity": gamma_zero,
        "non_gamma_zero_multiplicity": nongamma_zero,
        "bloch_zero_total": bloch_total,
        "direct_zero_count": direct_zero_count,
        "zero_eigenvalue_residual_below_tolerance": (
            max(direct_zero_residual, bloch_residual) <= tolerance
        ),
        "maximum_annihilation_residual": max(cls_residual, line_residual, global_residual),
        "checks": checks,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "L",
        "cells",
        "sites",
        "incidence_rank",
        "face_cls_count",
        "face_cls_span_rank",
        "face_cls_dependency_count",
        "noncontractible_line_count",
        "cls_plus_lines_rank",
        "edge_center_zero_dimension",
        "global_vertex_state_count",
        "vertex_zero_dimension",
        "complete_constructed_rank",
        "full_zero_dimension",
        "nominal_flat_band_states",
        "gamma_touching_extra_zero_modes",
        "gamma_zero_multiplicity",
        "non_gamma_zero_multiplicity",
        "zero_eigenvalue_residual_below_tolerance",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze(graph_path: Path, sizes: list[int], tolerance: float) -> dict[str, Any]:
    graph_path = graph_path.resolve()
    graph = read_graph(graph_path)
    rows = [analyze_size(graph, size, tolerance) for size in sizes]
    return {
        "schema_version": 1,
        "status": "verified_single_particle_graph_result",
        "scope": (
            "均匀最近邻紧束缚邻接矩阵的有限周期零能子空间；"
            "不外推到零场相互作用 Heisenberg 动力学谱"
        ),
        "graph": graph_path.relative_to(ROOT).as_posix(),
        "graph_sha256": sha256_file(graph_path),
        "boundary_condition": "L x L periodic torus",
        "numerical_zero_tolerance": tolerance,
        "rank_certificate": {
            "field": f"GF({RANK_PRIME})",
            "logic": (
                "有限域非零子式给出有理数域秩下界；二分关联矩阵、CLS 依赖关系"
                "与零空间维数给出同值上界，因此表中秩为精确整数秩"
            ),
        },
        "analytic_counting": {
            "cells": "M=L^2",
            "incidence_rank": "2M-1",
            "face_cls_span_rank": "M-1",
            "edge_center_zero_dimension": "M+1=(M-1)+2",
            "vertex_zero_dimension": "1",
            "full_zero_dimension": "M+2",
            "gamma_touching_extra_zero_modes": "(M+2)-M=2",
            "decomposition_beyond_cls": (
                "两个非收缩线态补全边心零空间；一个顶点全局态补全总零空间"
            ),
        },
        "state_construction": {
            "face_cls": (
                "每个蜂窝六边形的六个边心点依次取 +1,-1,+1,-1,+1,-1；"
                "所有 M 个面态之和为零"
            ),
            "line_a2": "固定 x，3 号边心取 +1、4 号边心取 -1，沿 y 绕环面一周",
            "line_a1": "固定 y，4 号边心取 +1、5 号边心取 -1，沿 x 绕环面一周",
            "global_vertex": "所有 1 号顶点取 +1，所有 2 号顶点取 -1，边心振幅为零",
        },
        "sizes": rows,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--sizes", nargs="+", type=int, default=[2, 3, 4, 5, 6, 8])
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.sizes or len(set(args.sizes)) != len(args.sizes) or any(size < 2 for size in args.sizes):
        raise SystemExit("尺寸必须是不重复且至少为 2 的整数")
    if args.tolerance <= 0.0:
        raise SystemExit("容差必须为正")
    report = analyze(args.graph, args.sizes, args.tolerance)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "summary.json"
    csv_path = args.output_dir / "size_counts.csv"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_csv(csv_path, report["sizes"])
    for row in report["sizes"]:
        print(
            f"L={row['L']}: CLS 秩={row['face_cls_span_rank']}, "
            f"边心零模={row['edge_center_zero_dimension']}, 总零模={row['full_zero_dimension']}"
        )
    print(f"结果：{json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
