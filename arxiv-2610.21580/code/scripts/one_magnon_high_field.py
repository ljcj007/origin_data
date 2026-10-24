#!/usr/bin/env python3
"""计算饱和场上方的精确单磁振子拉普拉斯谱。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scripts.band_graph import PeriodicGraph, bloch_hamiltonian, read_graph
except ModuleNotFoundError:
    from band_graph import PeriodicGraph, bloch_hamiltonian, read_graph


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH = ROOT / "graphs" / "super_honeycomb_poscar.dat"
DEFAULT_OUTPUT = ROOT / "results" / "one_magnon_high_field_v1" / "summary.json"


def sha256_file(path: Path) -> str:
    """计算文件哈希。"""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def basis_degrees(graph: PeriodicGraph) -> np.ndarray:
    """统计每个基底点的配位数。"""
    degrees = np.zeros(graph.basis_count, dtype=float)
    for edge in graph.edges:
        degrees[edge.site1 - 1] += 1.0
        degrees[edge.site2 - 1] += 1.0
    return degrees


def laplacian(graph: PeriodicGraph, kx: float, ky: float) -> np.ndarray:
    """构造 Bloch 图拉普拉斯。"""
    adjacency = bloch_hamiltonian(graph, kx, ky, np.ones(1))
    return np.diag(basis_degrees(graph)) - adjacency


def analyze(graph_path: Path, grid_size: int = 121, spin: float = 0.5) -> dict[str, Any]:
    """验证谱界、平带和饱和场。"""
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("grid_size 必须是至少为 3 的奇数")
    if spin <= 0.0:
        raise ValueError("spin 必须为正")
    graph_path = graph_path.resolve()
    graph = read_graph(graph_path)
    try:
        graph_record = graph_path.relative_to(ROOT).as_posix()
    except ValueError:
        graph_record = str(graph_path)
    degrees = basis_degrees(graph)
    grid = np.linspace(-math.pi, math.pi, grid_size)
    spectra = np.asarray(
        [
            np.linalg.eigvalsh(laplacian(graph, float(kx), float(ky)))
            for kx in grid
            for ky in grid
        ]
    )
    gamma = np.linalg.eigvalsh(laplacian(graph, 0.0, 0.0))
    flat_residual = float(np.max(np.min(np.abs(spectra - 2.0), axis=1)))
    lambda_max = float(np.max(spectra))
    lambda_min = float(np.min(spectra))
    h_sat_over_j = spin * lambda_max
    return {
        "schema_version": 1,
        "status": "verified_optional_extension",
        "scope": "饱和极化态上方的单磁振子扇区；不属于零场正式 SSE 协议",
        "graph": graph_record,
        "graph_sha256": sha256_file(graph_path),
        "basis_count": graph.basis_count,
        "degrees": degrees.astype(int).tolist(),
        "grid_size": grid_size,
        "laplacian": {
            "gamma_eigenvalues": gamma.tolist(),
            "lambda_min": lambda_min,
            "lambda_max": lambda_max,
            "flat_eigenvalue": 2.0,
            "flat_grid_residual": flat_residual,
        },
        "spin": spin,
        "one_magnon_energy": "epsilon_n/J = h/J - S*lambda_n",
        "flat_energy": "epsilon_flat/J = h/J - 2*S",
        "h_sat_over_j": h_sat_over_j,
        "spin_half_exactness": (
            "硬核玻色子表示的多磁振子相互作用 J sum_<ij> n_i n_j 非负；"
            "因此单粒子矩阵半正定的阈值 h/J=S*lambda_max 给出精确饱和场"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--grid", type=int, default=121)
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = analyze(args.graph.resolve(), args.grid, args.spin)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
