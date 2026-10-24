"""强星团解耦极限的热力学参照，不运行晶格 SSE。"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from check_super_honeycomb_extensions import spin_operators


LEVELS = np.array([-1.25, -0.75, 0.25, 0.75])
DEGENERACIES = np.array([3, 2, 6, 5])
TOTAL_SPINS = np.array([1, 0, 1, 2])
EDGE_TOTAL_Z = np.array([5/18, 0, 1/9, 1/2])


def reference(temperature, ja=1.0):
    if not math.isfinite(temperature) or not math.isfinite(ja) or temperature <= 0 or ja <= 0:
        raise ValueError('温度和强交换必须为有限正数')
    energies = ja * LEVELS
    excitation = energies - energies.min()
    weights = DEGENERACIES * np.exp(-excitation / temperature)
    partition = float(weights.sum())
    probabilities = weights / partition
    mean_excitation = float(probabilities @ excitation)
    energy = float(probabilities @ energies)
    variance = float(probabilities @ (excitation - mean_excitation)**2)
    q = float(probabilities @ EDGE_TOTAL_Z)
    return dict(T=temperature, JA=ja, energy_per_site=energy/5,
                specific_heat_per_site=variance/(5*temperature**2),
                entropy_per_site=(math.log(2*partition)+mean_excitation/temperature)/5,
                T_chi=(float(probabilities @ (TOTAL_SPINS*(TOTAL_SPINS+1)/3))+0.25)/5,
                strong_bond_correlation=energy/(3*ja),
                excited_star_probability=float(probabilities[1:].sum()),
                edge_total_z=q, d_T_chi_d_JB_at_zero=-3*q/(10*temperature))


def matrices():
    star = spin_operators()
    hstar = sum(star[0][a] @ star[j][a] for j in (1, 2, 3) for a in range(3))
    ops = [[np.kron(x, np.eye(2)) for x in site] for site in star]
    free = [np.kron(np.eye(16), x) for x in
            [np.array([[0, 1], [1, 0]])/2,
             np.array([[0, -1j], [1j, 0]])/2, np.diag([1, -1])/2]]
    mz = sum(site[2] for site in ops) + free[2]
    single_weak_bond = sum(ops[1][a] @ free[a] for a in range(3))
    return np.kron(hstar, np.eye(2)), mz @ mz, single_weak_bond


def thermal_trace(h, mz2, temperature):
    if h.shape != (32, 32) or np.max(np.abs(h-h.conj().T)) > 1e-13:
        raise ValueError('小矩阵维数或厄米性错误')
    values, vectors = np.linalg.eigh(h)
    residual = float(np.max(np.linalg.norm(h@vectors-vectors*values, axis=0)))
    if residual > 1e-12:
        raise ValueError('本征对残差超限')
    excitation = values-values.min()
    weights = np.exp(-excitation/temperature)
    probabilities = weights/weights.sum()
    mean_excitation = float(probabilities@excitation)
    diagonal = np.diag(vectors.conj().T@mz2@vectors).real
    return dict(energy_per_site=float(probabilities@values)/5,
                specific_heat_per_site=float(probabilities@(excitation-mean_excitation)**2)/(5*temperature**2),
                entropy_per_site=(math.log(float(weights.sum()))+mean_excitation/temperature)/5,
                T_chi=float(probabilities@diagonal)/5, eigenpair_residual=residual)


def verify():
    h, mz2, weak = matrices()
    temperatures = [0.01, 0.03, 0.1, 0.3, 1., 3., 10.]
    checks = []
    for temperature in temperatures:
        analytical = reference(temperature)
        direct = thermal_trace(h, mz2, temperature)
        error = max(abs(direct[k]-analytical[k]) for k in direct if k != 'eigenpair_residual')
        step = 1e-5*temperature
        # 三条等价弱键的一阶导数可相加；不近似高阶连通项。
        derivative = 3*(thermal_trace(h+step*weak, mz2, temperature)['T_chi']-
                        thermal_trace(h-step*weak, mz2, temperature)['T_chi'])/(2*step)
        derivative_error = abs(derivative-analytical['d_T_chi_d_JB_at_zero'])
        if error > 1e-11 or derivative_error > 3e-7:
            raise ValueError('解析式与独立热迹不一致')
        checks.append(dict(T=temperature, thermal_error=error,
                           weak_derivative_error=derivative_error,
                           eigenpair_residual=direct['eigenpair_residual']))
    return dict(status='passed', dimension=32, all_magnetization_sectors=True,
                weak_derivative_scope='single-bond derivative times three at JB=0, not finite-JB lattice ED',
                checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    validation = verify()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    rows = [reference(t) for t in [0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1., 3., 10.]]
    with (args.output_dir/'decoupled_reference.csv').open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    validation['provenance'] = {str(p.name): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__), Path(__file__).with_name('check_super_honeycomb_extensions.py')]}
    validation['scope'] = 'exact JB=0 thermodynamics and derivative; not new SSE data'
    with (args.output_dir/'validation.json').open('x', encoding='utf-8') as handle:
        json.dump(validation, handle, ensure_ascii=False, indent=2, allow_nan=False)
    print(json.dumps(validation, ensure_ascii=False))


if __name__ == '__main__':
    main()
