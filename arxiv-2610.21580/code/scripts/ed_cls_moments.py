"""检验 CLS 自旋激发的精确能量矩，不做谱线展宽。"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

import ed_graph_heisenberg as base


def sector(nsite, bonds, nup):
    diagonal, rows, columns, states = base.make_hamiltonian(nsite, bonds, nup)
    assert len(states) == math.comb(nsite, nup)
    matrix = coo_matrix((np.full(len(rows), 0.5), (rows, columns)),
                        shape=(len(states), len(states))).tocsr() + diags(diagonal)
    difference = matrix - matrix.T
    assert difference.nnz == 0 or np.max(np.abs(difference.data)) < 1e-14
    return matrix.tocsr(), states


def apply_spin(states, vector, target_states, weights, channel):
    output = np.zeros(len(target_states), dtype=np.result_type(vector, weights))
    index = {int(state): i for i, state in enumerate(target_states)}
    for site, weight in enumerate(weights):
        if weight == 0:
            continue
        up = ((states >> np.uint64(site)) & np.uint64(1)).astype(bool)
        if channel == "z":
            assert np.array_equal(states, target_states)
            output += weight * (up.astype(float) - 0.5) * vector
        else:
            chosen = np.flatnonzero(up if channel == "minus" else ~up)
            flipped = states[chosen] ^ np.uint64(1 << site)
            rows = np.fromiter((index[int(s)] for s in flipped), dtype=np.int64)
            np.add.at(output, rows, weight * vector[chosen])
    return output


def cls_weights(lx, ly, x=0, y=0):
    weights = np.zeros(5 * lx * ly)
    for local, dx, dy, sign in [(2,0,0,1), (3,0,0,-1), (4,-1,0,1),
                               (2,-1,0,-1), (3,-1,-1,1), (4,-1,-1,-1)]:
        weights[5 * (((x+dx) % lx) + lx * ((y+dy) % ly)) + local] += sign
    return weights


def moments(matrix, vector, energy0, ground=None):
    raw_weight = float(np.vdot(vector, vector).real)
    elastic_weight = 0.0
    if ground is not None:
        overlap = np.vdot(ground, vector)
        elastic_weight = float(abs(overlap)**2)
        vector = vector - overlap * ground
    weight = float(np.vdot(vector, vector).real)
    result = dict(raw_weight=raw_weight, elastic_weight=elastic_weight,
                  inelastic_weight=weight)
    if weight < 1e-20:
        result["null_inelastic_state"] = True
        return result
    vector = vector / np.sqrt(weight)
    shifted = matrix @ vector - energy0 * vector
    mean = float(np.vdot(vector, shifted).real)
    variance = float(np.vdot(shifted - mean * vector, shifted - mean * vector).real)
    result.update(mean_energy=mean, second_moment=float(np.vdot(shifted, shifted).real),
                  energy_variance=variance, energy_std=np.sqrt(variance),
                  null_inelastic_state=False)
    return result


def calculate(graph, lx, ly):
    nsite, bonds, eta, nplus, nminus = base.make_lattice(graph, lx, ly)
    if nsite > 20 or nsite < 10:
        raise ValueError("本检验仅允许 10<=N<=20")
    spin = abs(nplus - nminus) / 2
    nup = int(nsite / 2 + spin)
    matrix, states = sector(nsite, bonds, nup)
    seed = 20260905 + nup
    print(f"N={nsite}, Sz={spin}, D={len(states)}, Hermiticity=passed", flush=True)
    values, vectors = eigsh(matrix, k=1, which="SA", tol=1e-12, maxiter=1200,
                           v0=np.random.default_rng(seed).normal(size=len(states)))
    energy0, ground = float(values[0]), vectors[:, 0]
    residual = float(np.linalg.norm(matrix @ ground - energy0 * ground))
    assert residual < 1e-9
    print(f"E0={energy0:.14f}, residual={residual:.3e}", flush=True)
    adjacency = np.zeros((nsite, nsite))
    for i, j in bonds:
        adjacency[i, j] = adjacency[j, i] = 1
    weights = cls_weights(lx, ly)
    assert np.linalg.norm(weights) > 0 and np.linalg.norm(adjacency @ weights) == 0
    weights /= np.linalg.norm(weights)
    probes = {"local_cls": weights}
    # 2x2 的离散动量只需实数相位。
    if lx == ly == 2:
        for kx, ky in [(0,0), (1,0), (0,1), (1,1)]:
            coefficients = sum(((-1)**(kx*x + ky*y)) * cls_weights(lx,ly,x,y)
                               for x in range(lx) for y in range(ly))
            norm = np.linalg.norm(coefficients)
            if norm > 0:
                probes[f"cls_q_{kx}pi_{ky}pi"] = coefficients / norm
            else:
                assert kx == ky == 0
    channels, controls = {}, {}
    for channel, delta in [("minus",-1), ("z",0), ("plus",1)]:
        target_h, target = (matrix, states) if delta == 0 else sector(nsite,bonds,nup+delta)
        uniform = apply_spin(states,ground,target,np.ones(nsite),channel)
        control = moments(target_h,uniform,energy0)
        controls[channel] = control
        if channel == "minus":
            assert abs(control["raw_weight"] - 2*spin) < 1e-9
            multiplet = uniform / np.linalg.norm(uniform)
            assert np.linalg.norm(target_h @ multiplet - energy0*multiplet) < 1e-9
        elif channel == "z":
            assert abs(control["raw_weight"] - spin**2) < 1e-9
            multiplet = ground
        else:
            assert control["raw_weight"] < 1e-18
            multiplet = None
        channels[channel] = {name: moments(target_h,
            apply_spin(states,ground,target,w,channel),energy0,multiplet)
            for name,w in probes.items()}
        print(f"channel={channel}, D={len(target)}, {channels[channel]['local_cls']}",flush=True)
    return dict(schema_version=1,status="completed", lx=lx,ly=ly,nsite=nsite,nbond=len(bonds),
        hamiltonian="H=sum_<ij> S_i dot S_j; S=1/2; J=1; h=0; PBC",
        ground_spin=spin, ground_magnetization=spin, ground_energy=energy0,
        ground_residual=residual, solver=dict(method="eigsh",tol=1e-12,maxiter=1200,seed=seed),
        numpy_version=np.__version__,scipy_version=scipy.__version__,
        graph_sha256=base.sha256_file(graph),source_sha256=base.sha256_file(Path(__file__)),
        base_source_sha256=base.sha256_file(Path(base.__file__)),
        cls_weights=weights.tolist(),adjacency_cls_residual=float(np.linalg.norm(adjacency@weights)),
        controls=controls,channels=channels,
        interpretation="Spectral moments only; nonzero variance excludes a single-energy state for this probe, not all flat branches.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph",type=Path,required=True)
    parser.add_argument("--lx",type=int,default=2)
    parser.add_argument("--ly",type=int,default=2)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    base.require_new_output(args.output)
    result=calculate(args.graph,args.lx,args.ly)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x",encoding="utf-8") as handle:
        json.dump(result,handle,ensure_ascii=False,indent=2,allow_nan=False)
        handle.write("\n")
