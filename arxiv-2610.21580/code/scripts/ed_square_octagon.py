"""N=16 square-octagon 海森堡模型的精确对角化校验。"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np


def site(cell: int, local_site: int) -> int:
    return 4 * cell + local_site


def make_lattice(lx: int, ly: int) -> tuple[int, list[tuple[int, int]], np.ndarray]:
    ncell = lx * ly
    nsite = 4 * ncell
    bonds: list[tuple[int, int]] = []
    eta = np.empty(nsite, dtype=np.int8)

    def cell(x: int, y: int) -> int:
        return (x % lx) + (y % ly) * lx

    for y in range(ly):
        for x in range(lx):
            c = cell(x, y)
            cr = cell(x + 1, y)
            cu = cell(x, y + 1)
            bonds.extend(
                [
                    (site(c, 0), site(c, 1)),
                    (site(c, 1), site(c, 2)),
                    (site(c, 2), site(c, 3)),
                    (site(c, 3), site(c, 0)),
                    (site(c, 1), site(cr, 3)),
                    (site(c, 2), site(cu, 0)),
                ]
            )
            for a in range(4):
                eta[site(c, a)] = 1 if (a + 1 + x + y) % 2 == 0 else -1

    if len(set(tuple(sorted(b)) for b in bonds)) != len(bonds):
        raise ValueError("存在重复键")
    return nsite, bonds, eta


def make_basis(nsite: int, nup: int) -> tuple[np.ndarray, dict[int, int]]:
    states = np.fromiter(
        (sum(1 << i for i in ups) for ups in combinations(range(nsite), nup)),
        dtype=np.uint64,
    )
    index = {int(state): i for i, state in enumerate(states)}
    return states, index


def make_hamiltonian(
    nsite: int, bonds: list[tuple[int, int]], nup: int
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    states, index = make_basis(nsite, nup)
    diagonal = np.zeros(len(states))
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []

    for i, j in bonds:
        up_i = (states >> np.uint64(i)) & np.uint64(1)
        up_j = (states >> np.uint64(j)) & np.uint64(1)
        opposite = up_i != up_j
        diagonal += np.where(opposite, -0.25, 0.25)
        source = np.flatnonzero(opposite)
        flipped = states[source] ^ np.uint64(1 << i) ^ np.uint64(1 << j)
        target = np.fromiter((index[int(s)] for s in flipped), dtype=np.int64)
        rows.append(target)
        cols.append(source)

    return (diagonal, np.concatenate(rows), np.concatenate(cols)), states


def matvec(matrix, vector: np.ndarray) -> np.ndarray:
    diagonal, rows, cols = matrix
    result = diagonal * vector
    np.add.at(result, rows, 0.5 * vector[cols])
    return result


def lowest_state(matrix, maxiter: int = 160) -> tuple[float, np.ndarray, float]:
    """完全重正交化 Lanczos，只求最低本征态。"""
    dimension = len(matrix[0])
    rng = np.random.default_rng(20260902)
    q = rng.normal(size=dimension)
    q /= np.linalg.norm(q)
    basis = np.empty((dimension, maxiter))
    alphas: list[float] = []
    betas: list[float] = []

    for k in range(maxiter):
        basis[:, k] = q
        z = matvec(matrix, q)
        if k > 0:
            z -= betas[k - 1] * basis[:, k - 1]
        alpha = float(np.dot(q, z))
        z -= alpha * q

        block = basis[:, : k + 1]
        z -= block @ (block.T @ z)
        z -= block @ (block.T @ z)
        beta = float(np.linalg.norm(z))
        alphas.append(alpha)
        if beta < 1.0e-13:
            break
        betas.append(beta)
        q = z / beta

    m = len(alphas)
    tridiagonal = np.diag(alphas)
    if m > 1:
        offdiag = np.asarray(betas[: m - 1])
        tridiagonal += np.diag(offdiag, 1) + np.diag(offdiag, -1)
    values, vectors = np.linalg.eigh(tridiagonal)
    coeff = vectors[:, 0]
    state = basis[:, :m] @ coeff
    residual = 0.0 if len(betas) < m else abs(betas[m - 1] * coeff[-1])
    return float(values[0]), state, float(residual)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lx", type=int, default=2)
    parser.add_argument("--ly", type=int, default=2)
    parser.add_argument("--output", type=Path, default=Path("results/ed_L2.json"))
    args = parser.parse_args()

    nsite, bonds, eta = make_lattice(args.lx, args.ly)
    if nsite > 20:
        raise ValueError("教学版 ED 仅用于 N<=20 的校验")

    h0, states0 = make_hamiltonian(nsite, bonds, nsite // 2)
    e0, psi0, residual0 = lowest_state(h0)
    h1, _ = make_hamiltonian(nsite, bonds, nsite // 2 + 1)
    e1, _, residual1 = lowest_state(h1)

    ms_values = np.empty(len(states0))
    for k, raw_state in enumerate(states0):
        state = int(raw_state)
        spins = np.array([1 if (state >> i) & 1 else -1 for i in range(nsite)])
        ms_values[k] = 0.5 * np.dot(eta, spins)
    staggered_m2 = 3.0 * np.dot(np.abs(psi0) ** 2, ms_values**2) / nsite**2

    result = {
        "lx": args.lx,
        "ly": args.ly,
        "nsite": nsite,
        "nbond": len(bonds),
        "ground_state_energy": e0,
        "ground_state_energy_per_site": e0 / nsite,
        "triplet_gap": e1 - e0,
        "ground_state_staggered_m2": staggered_m2,
        "lanczos_residual_sz0": residual0,
        "lanczos_residual_sz1": residual1,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
