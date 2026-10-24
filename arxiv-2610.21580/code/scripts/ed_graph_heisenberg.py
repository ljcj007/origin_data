#!/usr/bin/env python3
"""对周期格点图的自旋 1/2 海森堡模型做小尺寸 ED。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def require_new_output(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"输出已存在，拒绝覆盖：{path}")


def read_graph(path: Path) -> tuple[int, list[tuple[int, int, int, int]], np.ndarray, np.ndarray]:
    lines = [
        content.split()
        for raw in path.read_text(encoding="utf-8").splitlines()
        if (content := raw.split("#", 1)[0].strip())
    ]
    if not lines or len(lines[0]) != 3:
        raise ValueError("格点图首行格式错误")
    n_basis, n_edges, dimension = (int(value) for value in lines[0][:3])
    if n_basis < 1 or n_edges < 1 or dimension != 3:
        raise ValueError("格点图头信息错误")
    if len(lines) != n_edges + 3:
        raise ValueError("格点图记录数量错误")
    edges = []
    for tokens in lines[1 : 1 + n_edges]:
        if len(tokens) != 6:
            raise ValueError("格点图键记录格式错误")
        site1, site2, dx, dy, dz, bond_type = (int(value) for value in tokens)
        if not 1 <= site1 <= n_basis or not 1 <= site2 <= n_basis:
            raise ValueError("格点图含越界编号")
        if dz != 0:
            raise ValueError("当前脚本仅支持 xy 平面内的键")
        if bond_type != 1:
            raise ValueError("教学版 ED 仅支持一种键")
        edges.append((site1 - 1, site2 - 1, dx, dy))
    eta_tokens = lines[1 + n_edges]
    flip_tokens = lines[2 + n_edges]
    if len(eta_tokens) != n_basis or len(flip_tokens) != 3:
        raise ValueError("格点图子格记录格式错误")
    eta = np.asarray([int(value) for value in eta_tokens], dtype=np.int8)
    flips = np.asarray([int(value) for value in flip_tokens], dtype=np.int8)
    if len(eta) != n_basis or not np.all(np.abs(eta) == 1):
        raise ValueError("格点图没有有效二分子格标记")
    if not np.all((flips == 0) | (flips == 1)):
        raise ValueError("子格平移相位必须为 0 或 1")
    return n_basis, edges, eta, flips[:2]


def make_lattice(
    graph: Path, lx: int, ly: int
) -> tuple[int, list[tuple[int, int]], np.ndarray, int, int]:
    if lx < 1 or ly < 1:
        raise ValueError("Lx 和 Ly 必须为正整数")
    n_basis, edge_defs, eta_basis, flips = read_graph(graph)
    nsite = n_basis * lx * ly
    bonds: list[tuple[int, int]] = []
    eta = np.empty(nsite, dtype=np.int8)

    def site(x: int, y: int, local: int) -> int:
        return n_basis * ((x % lx) + (y % ly) * lx) + local

    for y in range(ly):
        for x in range(lx):
            parity = (x * int(flips[0]) + y * int(flips[1])) % 2
            for local in range(n_basis):
                eta[site(x, y, local)] = eta_basis[local] * (-1 if parity else 1)
            for local1, local2, dx, dy in edge_defs:
                i, j = site(x, y, local1), site(x + dx, y + dy, local2)
                if i == j:
                    raise ValueError("周期化后出现自键")
                bonds.append((i, j))
    undirected = {(min(i, j), max(i, j)) for i, j in bonds}
    if len(undirected) != len(bonds):
        raise ValueError("周期化后出现重复无向键")
    if any(eta[i] * eta[j] != -1 for i, j in bonds):
        raise ValueError("有限晶格不是二分图")
    n_plus = int(np.count_nonzero(eta == 1))
    n_minus = nsite - n_plus
    return nsite, bonds, eta, n_plus, n_minus


def make_basis(nsite: int, nup: int) -> tuple[np.ndarray, dict[int, int]]:
    if not 0 <= nup <= nsite:
        raise ValueError("nup 超出允许范围")
    states = np.fromiter(
        (sum(1 << site for site in ups) for ups in combinations(range(nsite), nup)),
        dtype=np.uint64,
    )
    return states, {int(state): index for index, state in enumerate(states)}


def make_hamiltonian(
    nsite: int, bonds: list[tuple[int, int]], nup: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    states, index = make_basis(nsite, nup)
    diagonal = np.zeros(len(states), dtype=np.float64)
    row_blocks: list[np.ndarray] = []
    column_blocks: list[np.ndarray] = []
    for i, j in bonds:
        up_i = (states >> np.uint64(i)) & np.uint64(1)
        up_j = (states >> np.uint64(j)) & np.uint64(1)
        opposite = up_i != up_j
        diagonal += np.where(opposite, -0.25, 0.25)
        columns = np.flatnonzero(opposite)
        flipped = states[columns] ^ np.uint64(1 << i) ^ np.uint64(1 << j)
        rows = np.fromiter((index[int(state)] for state in flipped), dtype=np.int64)
        row_blocks.append(rows)
        column_blocks.append(columns)
    return diagonal, np.concatenate(row_blocks), np.concatenate(column_blocks), states


def lowest_state(
    nsite: int,
    bonds: list[tuple[int, int]],
    nup: int,
    tol: float,
    maxiter: int,
    return_vector: bool,
) -> tuple[float, np.ndarray | None, np.ndarray, float, float]:
    diagonal, rows, columns, states = make_hamiltonian(nsite, bonds, nup)

    def matvec(vector: np.ndarray) -> np.ndarray:
        result = diagonal * vector
        np.add.at(result, rows, 0.5 * vector[columns])
        return result

    operator = LinearOperator((len(states), len(states)), matvec=matvec, dtype=np.float64)
    rng = np.random.default_rng(20260903 + nup)
    values, vectors = eigsh(
        operator,
        k=1,
        which="SA",
        v0=rng.normal(size=len(states)),
        tol=tol,
        maxiter=maxiter,
    )
    vector = vectors[:, 0]
    residual = float(np.linalg.norm(matvec(vector) - values[0] * vector))
    state_norm = float(np.dot(vector, vector))
    return float(values[0]), vector if return_vector else None, states, residual, state_norm


def ladder_norm(
    states: np.ndarray,
    vector: np.ndarray,
    weights: np.ndarray,
    nsite: int,
    raise_spin: bool,
) -> float:
    nup = bin(int(states[0])).count("1")
    target_nup = nup + (1 if raise_spin else -1)
    target_states, target_index = make_basis(nsite, target_nup)
    output = np.zeros(len(target_states), dtype=np.float64)
    for site in range(nsite):
        occupied = ((states >> np.uint64(site)) & np.uint64(1)).astype(bool)
        selected = np.flatnonzero(~occupied if raise_spin else occupied)
        changed = states[selected] ^ np.uint64(1 << site)
        targets = np.fromiter((target_index[int(state)] for state in changed), dtype=np.int64)
        np.add.at(output, targets, float(weights[site]) * vector[selected])
    return float(np.dot(output, output))


def vector_square(
    states: np.ndarray, vector: np.ndarray, weights: np.ndarray, nsite: int
) -> float:
    mz = np.zeros(len(states), dtype=np.float64)
    for site in range(nsite):
        occupied = ((states >> np.uint64(site)) & np.uint64(1)).astype(np.float64)
        mz += float(weights[site]) * (occupied - 0.5)
    z2 = float(np.dot(vector * vector, mz * mz))
    plus = ladder_norm(states, vector, weights, nsite, True)
    minus = ladder_norm(states, vector, weights, nsite, False)
    return z2 + 0.5 * (plus + minus)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=Path("graphs/super_honeycomb_poscar.dat"))
    parser.add_argument("--lx", type=int, default=2)
    parser.add_argument("--ly", type=int, default=2)
    parser.add_argument("--tol", type=float, default=1.0e-11)
    parser.add_argument("--maxiter", type=int, default=800)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/ed_super_honeycomb_v1/result.json"),
    )
    args = parser.parse_args()

    if not 0.0 < args.tol <= 1.0e-9:
        raise ValueError("tol 必须位于 (0, 1e-9]")
    if args.maxiter < 1:
        raise ValueError("maxiter 必须为正整数")
    require_new_output(args.output)

    nsite, bonds, eta, n_plus, n_minus = make_lattice(args.graph, args.lx, args.ly)
    if nsite > 20:
        raise ValueError("教学版 ED 仅用于 N<=20")
    if nsite % 2:
        raise ValueError("当前脚本要求偶数个格点")
    ground_spin = 0.5 * abs(n_plus - n_minus)
    max_sector = int(ground_spin) + 1
    sectors = []
    ground_vector = None
    ground_states = None
    for magnetization in range(max_sector + 1):
        nup = nsite // 2 + magnetization
        energy, vector, states, residual, state_norm = lowest_state(
            nsite,
            bonds,
            nup,
            args.tol,
            args.maxiter,
            magnetization == 0,
        )
        sectors.append(
            {
                "magnetization": magnetization,
                "nup": nup,
                "random_seed": 20260903 + nup,
                "dimension": len(states),
                "expected_dimension": math.comb(nsite, nup),
                "energy": energy,
                "residual": residual,
                "state_norm": state_norm,
            }
        )
        if magnetization == 0:
            ground_vector, ground_states = vector, states

    assert ground_vector is not None and ground_states is not None
    uniform_s2 = vector_square(ground_states, ground_vector, np.ones(nsite), nsite)
    staggered_s2 = vector_square(ground_states, ground_vector, eta.astype(float), nsite)
    energy0 = float(sectors[0]["energy"])
    result = {
        "schema_version": 2,
        "status": "completed",
        "graph": str(args.graph),
        "graph_sha256": sha256_file(args.graph),
        "source_sha256": sha256_file(Path(__file__)),
        "hamiltonian": "H = sum_<ij> S_i dot S_j, J=1",
        "solver": {
            "method": "scipy.sparse.linalg.eigsh",
            "tol": args.tol,
            "maxiter": args.maxiter,
            "random_seed_rule": "20260903+nup",
        },
        "lx": args.lx,
        "ly": args.ly,
        "nsite": nsite,
        "nbond": len(bonds),
        "n_plus": n_plus,
        "n_minus": n_minus,
        "lieb_mattis_ground_spin": ground_spin,
        "sector_energies": sectors,
        "ground_multiplet_energy_spread": max(
            abs(float(item["energy"]) - energy0)
            for item in sectors
            if int(item["magnetization"]) <= int(ground_spin)
        ),
        "gap_to_next_spin_sector": float(sectors[-1]["energy"]) - energy0,
        "ground_state_energy": energy0,
        "ground_state_energy_per_site": energy0 / nsite,
        "maximum_lanczos_residual": max(float(item["residual"]) for item in sectors),
        "uniform_spin_squared": uniform_s2,
        "lieb_mattis_spin_squared": ground_spin * (ground_spin + 1.0),
        "ferrimagnetic_m2": uniform_s2 / nsite**2,
        "staggered_m2": staggered_s2 / nsite**2,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
