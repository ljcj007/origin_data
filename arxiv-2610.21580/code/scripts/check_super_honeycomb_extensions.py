#!/usr/bin/env python3
"""核对两类交换的平带、饱和场与四自旋星团极限。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from band_graph import GraphEdge, PeriodicGraph, bloch_hamiltonian, read_graph


ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT / "graphs/super_honeycomb_poscar.dat"
OUT = ROOT / "results/super_honeycomb_literature_extension_v1"


def spin_operators():
    """四个自旋的完整 16 维空间。"""
    local = [np.array([[0, 1], [1, 0]]) / 2,
             np.array([[0, -1j], [1j, 0]]) / 2,
             np.diag([1, -1]) / 2]
    ops = []
    for site in range(4):
        components = []
        for component in local:
            matrix = np.ones((1, 1))
            for j in range(4):
                matrix = np.kron(matrix, component if site == j else np.eye(2))
            components.append(matrix)
        ops.append(components)
    return ops


def star_check():
    ops = spin_operators()
    h = sum(ops[0][a] @ ops[j][a] for j in (1, 2, 3) for a in range(3))
    assert h.shape == (16, 16)
    assert np.max(np.abs(h - h.conj().T)) < 1e-12
    energies, vectors = np.linalg.eigh(h)
    expected = np.sort([-1.25] * 3 + [-0.75] * 2 + [0.25] * 6 + [0.75] * 5)
    ground = vectors[:, np.abs(energies + 1.25) < 1e-10]
    total = [sum(ops[j][a] for j in range(4)) for a in range(3)]
    projected = [ground.conj().T @ op @ ground for op in total]
    residuals = []
    for site, coefficient in enumerate([-0.25, 5/12, 5/12, 5/12]):
        for a in range(3):
            actual = ground.conj().T @ ops[site][a] @ ground
            residuals.append(float(np.linalg.norm(actual - coefficient * projected[a])))
    spin2 = ground.conj().T @ sum(op @ op for op in total) @ ground
    result = {
        "hilbert_dimension": 16,
        "route": "全空间 Kronecker 矩阵与 NumPy eigh；小星团代数核对",
        "energies_over_JA": energies.tolist(),
        "spectrum_error": float(np.max(np.abs(energies - expected))),
        "ground_multiplicity": ground.shape[1],
        "ground_energy_over_JA": float(energies[0]),
        "gap_over_JA": float(energies[3] - energies[0]),
        "spin2_error": float(np.linalg.norm(spin2 - 2 * np.eye(3))),
        "projection_max_residual": max(residuals),
        "central_projection": -0.25,
        "each_edge_projection": 5/12,
        "effective_exchange_over_JB": 5/12,
    }
    assert result["spectrum_error"] < 1e-12
    assert result["spin2_error"] < 1e-12
    assert result["projection_max_residual"] < 1e-12
    return result


def weighted_checks():
    source = read_graph(GRAPH)
    assert source.basis_count == 5 and len(source.edges) == 6
    assert all(e.site1 in (1, 2) and e.site2 in (3, 4, 5) for e in source.edges)
    graph = PeriodicGraph(5, tuple(GraphEdge(e.site1, e.site2, e.shift, e.site1)
                                  for e in source.edges), source.eta, source.phase)
    grid = np.linspace(-np.pi, np.pi, 31)
    rows = []
    for a, b in [(1., 0.1), (1., 0.25), (1., 0.5), (1., 1.), (1., 2.)]:
        degree = np.zeros(5)
        for edge in graph.edges:
            weight = (a, b)[edge.bond_type - 1]
            degree[edge.site1 - 1] += weight
            degree[edge.site2 - 1] += weight
        max_lap = -np.inf
        tb_res = lap_res = sw_res = 0.
        for k1 in grid:
            for k2 in grid:
                adjacency = bloch_hamiltonian(graph, k1, k2, np.array([a, b]))
                laplacian = np.diag(degree) - adjacency
                lap_values = np.linalg.eigvalsh(laplacian)
                max_lap = max(max_lap, lap_values[-1])
                tb_res = max(tb_res, np.min(np.abs(np.linalg.eigvalsh(adjacency))))
                lap_res = max(lap_res, np.min(np.abs(lap_values - (a + b))))
                c = adjacency[:2, 2:]
                dynamical = 0.5 * np.block([[np.diag(degree[:2]), -c],
                                           [c.conj().T, -np.diag(degree[2:])]])
                sw_res = max(sw_res, np.min(np.abs(np.linalg.eigvals(dynamical) + (a+b)/2)))
        predicted_max = 2*(a+b) + np.sqrt(4*a*a - 7*a*b + 4*b*b)
        row = {"JA": a, "JB": b, "grid": [31, 31],
               "tb_zero_residual": float(tb_res),
               "lswt_flat_energy": (a+b)/2,
               "lswt_flat_residual": float(sw_res),
               "laplacian_flat_eigenvalue": a+b,
               "laplacian_flat_residual": float(lap_res),
               "hsat_spin_half": float(predicted_max/2),
               "hsat_grid_error": float(abs(max_lap-predicted_max)/2),
               "flat_excitation_at_hsat": float(predicted_max/2-(a+b)/2)}
        assert max(tb_res, sw_res, lap_res, abs(max_lap-predicted_max)) < 1e-10
        rows.append(row)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 1, "status": "verified_analytic_extension",
              "date": "2026-09-05", "spin": 0.5,
              "boundary": "周期体相；四自旋星团为开放三臂星形",
              "scope": "新增解析检验，不是正式 SSE 扫描或零场动力学证明",
              "graph_sha256": hashlib.sha256(GRAPH.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "star": star_check(), "weighted_spectra": weighted_checks()}
    (OUT / "analytic_checks.json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
