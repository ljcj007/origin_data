#!/usr/bin/env python3
"""精确枚举周期图上的经典反铁磁 Ising 模型。"""

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
class Graph:
    nbasis: int
    edges: tuple[tuple[int, int, int, int], ...]
    eta_basis: tuple[int, ...]
    phase: tuple[int, int]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_graph(path: Path) -> Graph:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError("格点图为空")
    nbasis, nedge, dimension = map(int, rows[0])
    if dimension != 3 or len(rows) != nedge + 3:
        raise ValueError("格点图格式错误")
    edges = []
    for row in rows[1 : nedge + 1]:
        i, j, dx, dy, dz, bond_type = map(int, row)
        if not (1 <= i <= nbasis and 1 <= j <= nbasis):
            raise ValueError("格点编号越界")
        if dz != 0 or bond_type < 1:
            raise ValueError("仅支持二维正耦合键类型")
        edges.append((i, j, dx, dy))
    eta_basis = tuple(map(int, rows[nedge + 1]))
    phase3 = tuple(map(int, rows[nedge + 2]))
    if len(eta_basis) != nbasis or any(abs(value) != 1 for value in eta_basis):
        raise ValueError("eta 定义错误")
    if len(phase3) != 3 or phase3[2] != 0 or any(value not in (0, 1) for value in phase3):
        raise ValueError("平移相位定义错误")
    return Graph(nbasis, tuple(edges), eta_basis, (phase3[0], phase3[1]))


def build_supercell(graph: Graph, lx: int, ly: int) -> tuple[np.ndarray, np.ndarray]:
    if lx < 1 or ly < 1:
        raise ValueError("超胞尺寸必须为正")

    def site_index(x: int, y: int, basis: int) -> int:
        return basis - 1 + graph.nbasis * ((x % lx) + (y % ly) * lx)

    bonds: list[tuple[int, int]] = []
    eta = np.empty(graph.nbasis * lx * ly, dtype=np.int8)
    for y in range(ly):
        for x in range(lx):
            for basis in range(1, graph.nbasis + 1):
                site = site_index(x, y, basis)
                eta[site] = graph.eta_basis[basis - 1] * (-1) ** (
                    graph.phase[0] * x + graph.phase[1] * y
                )
            for i, j, dx, dy in graph.edges:
                bonds.append((site_index(x, y, i), site_index(x + dx, y + dy, j)))
    bond_array = np.asarray(bonds, dtype=np.int32)
    if np.any(bond_array[:, 0] == bond_array[:, 1]):
        raise ValueError("有限超胞产生自连接")
    if np.any(eta[bond_array[:, 0]] * eta[bond_array[:, 1]] != -1):
        raise ValueError("有限超胞不是二分图")
    canonical = {tuple(sorted(map(int, bond))) for bond in bond_array}
    if len(canonical) != len(bonds):
        raise ValueError("有限超胞产生重复键")
    return bond_array, eta


def enumerate_density(
    bonds: np.ndarray, eta: np.ndarray, chunk_size: int = 131072
) -> tuple[np.ndarray, np.ndarray]:
    nsite = len(eta)
    if nsite > 24:
        raise ValueError("精确枚举限制为 N<=24")
    if chunk_size < 1:
        raise ValueError("chunk_size 必须为正")
    nconfig = 1 << nsite
    chunks: list[np.ndarray] = []
    bit_positions = np.arange(nsite, dtype=np.uint32)
    for start in range(0, nconfig, chunk_size):
        states = np.arange(start, min(start + chunk_size, nconfig), dtype=np.uint32)
        spins = (((states[:, None] >> bit_positions) & 1) * 2 - 1).astype(np.int8)
        energy = np.zeros(len(states), dtype=np.int16)
        for site1, site2 in bonds:
            energy += spins[:, site1] * spins[:, site2]
        magnetization = spins.sum(axis=1, dtype=np.int16)
        staggered = (spins * eta).sum(axis=1, dtype=np.int16)
        chunks.append(np.column_stack((energy, magnetization, staggered)))
    states = np.concatenate(chunks, axis=0)
    density, degeneracy = np.unique(states, axis=0, return_counts=True)
    if int(degeneracy.sum()) != nconfig:
        raise RuntimeError("态数守恒检查失败")
    return density.astype(np.int16, copy=False), degeneracy.astype(np.int64, copy=False)


def thermal_observables(
    density: np.ndarray,
    degeneracy: np.ndarray,
    nsite: int,
    temperature: float,
    coupling: float = 1.0,
) -> dict[str, float]:
    if temperature <= 0.0 or coupling <= 0.0:
        raise ValueError("T 和 J 必须为正")
    energy = coupling * density[:, 0].astype(np.float64)
    magnetization = density[:, 1].astype(np.float64)
    staggered = density[:, 2].astype(np.float64)
    minimum = float(np.min(energy))
    weights = degeneracy.astype(np.float64) * np.exp(-(energy - minimum) / temperature)
    scaled_partition = float(np.sum(weights))

    def mean(values: np.ndarray) -> float:
        return float(np.dot(weights, values) / scaled_partition)

    mean_energy = mean(energy)
    mean_energy2 = mean(energy**2)
    mean_m = mean(magnetization)
    mean_m2 = mean(magnetization**2)
    mean_abs_ms = mean(np.abs(staggered))
    mean_ms2 = mean(staggered**2)
    mean_ms4 = mean(staggered**4)
    return {
        "energy_per_site": mean_energy / nsite,
        "specific_heat_per_site": (mean_energy2 - mean_energy**2) / (nsite * temperature**2),
        "abs_staggered_m": mean_abs_ms / nsite,
        "staggered_m2": mean_ms2 / nsite**2,
        "staggered_m4": mean_ms4 / nsite**4,
        "binder_staggered": 1.0 - mean_ms4 / (3.0 * mean_ms2**2),
        "uniform_susceptibility": (mean_m2 - mean_m**2) / (nsite * temperature),
        "connected_staggered_susceptibility": (mean_ms2 - mean_abs_ms**2)
        / (nsite * temperature),
        "mean_energy": mean_energy,
        "mean_energy_squared": mean_energy2,
        "mean_magnetization": mean_m,
        "mean_magnetization_squared": mean_m2,
        "mean_abs_staggered_magnetization": mean_abs_ms,
        "mean_staggered_magnetization_squared": mean_ms2,
        "mean_staggered_magnetization_fourth": mean_ms4,
        "log_partition_function": -minimum / temperature + math.log(scaled_partition),
    }


def write_density(path: Path, density: np.ndarray, degeneracy: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["energy_over_J", "magnetization", "staggered_magnetization", "degeneracy"])
        for values, count in zip(density, degeneracy, strict=True):
            writer.writerow([int(values[0]), int(values[1]), int(values[2]), int(count)])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, required=True)
    parser.add_argument("--lx", type=int, default=2)
    parser.add_argument("--ly", type=int, default=2)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--coupling", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    graph_path = args.graph.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    graph = read_graph(graph_path)
    bonds, eta = build_supercell(graph, args.lx, args.ly)
    density, degeneracy = enumerate_density(bonds, eta)
    nsite = len(eta)
    observables = thermal_observables(density, degeneracy, nsite, args.temperature, args.coupling)
    ground_mask = density[:, 0] == int(np.min(density[:, 0]))
    exact = {
        "schema_version": 1,
        "model": "classical_antiferromagnetic_ising",
        "conventions": {
            "hamiltonian": "H=+J*sum_<ij> s_i*s_j, J>0, h=0",
            "spin_values": [-1, 1],
            "staggered_magnetization": "Ms=sum_i eta_i*s_i",
            "normalization": {
                "energy_per_site": "<H>/N",
                "specific_heat_per_site": "(<H^2>-<H>^2)/(N*T^2)",
                "abs_staggered_m": "<|Ms|>/N",
                "staggered_m2": "<Ms^2>/N^2",
                "staggered_m4": "<Ms^4>/N^4",
                "binder_staggered": "1-<Ms^4>/(3*<Ms^2>^2)",
                "uniform_susceptibility": "(<M^2>-<M>^2)/(N*T)",
                "connected_staggered_susceptibility": "(<Ms^2>-<|Ms|>^2)/(N*T)",
            },
        },
        "lattice": {
            "graph": graph_path.as_posix(),
            "graph_sha256": sha256_file(graph_path),
            "boundary": "periodic_xy",
            "Lx": args.lx,
            "Ly": args.ly,
            "N": nsite,
            "Nb": len(bonds),
            "eta_plus": int(np.count_nonzero(eta == 1)),
            "eta_minus": int(np.count_nonzero(eta == -1)),
            "eta": [int(value) for value in eta],
            "bonds_zero_based": [[int(site1), int(site2)] for site1, site2 in bonds],
        },
        "parameters": {"temperature": args.temperature, "coupling": args.coupling, "field": 0.0},
        "enumeration": {
            "number_of_configurations": 1 << nsite,
            "density_rows": len(density),
            "minimum_energy_over_J": int(np.min(density[:, 0])),
            "ground_state_degeneracy": int(np.sum(degeneracy[ground_mask])),
            "global_spin_flip_symmetry": bool(
                abs(observables["mean_magnetization"]) < 1.0e-13
            ),
        },
        "observables": observables,
    }
    write_density(output_dir / "density_of_states.csv", density, degeneracy)
    (output_dir / "exact.json").write_text(
        json.dumps(exact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(output_dir / "exact.json")


if __name__ == "__main__":
    main()
