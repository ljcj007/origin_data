"""检验能量积分的熵估计与误差传播，不读取或生成晶格 QMC 数据。"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def entropy_weights(beta):
    beta = np.asarray(beta, dtype=float)
    if (beta.ndim != 1 or len(beta) < 2 or not np.all(np.isfinite(beta))
            or beta[0] <= 0 or np.any(np.diff(beta) <= 0)):
        raise ValueError('beta 必须是严格递增的有限正数，至少两点')
    integral = np.zeros((len(beta), len(beta)))
    for i, width in enumerate(np.diff(beta), start=1):
        integral[i] = integral[i-1]
        integral[i, i-1:i+1] += width/2
    return np.diag(beta)-integral


def reconstruct(beta, energy, tail_integral, covariance=None, tail_abs_bound=None):
    """按每原始自旋归一；统计协方差以固定高温尾部积分为条件。"""
    weights = entropy_weights(beta)
    energy = np.asarray(energy, dtype=float)
    if energy.shape != (len(weights),) or not np.all(np.isfinite(energy)):
        raise ValueError('能量必须与 beta 等长且全部有限')
    if not math.isfinite(tail_integral):
        raise ValueError('高温尾部积分必须有限')
    if tail_abs_bound is not None and (not math.isfinite(tail_abs_bound) or tail_abs_bound < 0):
        raise ValueError('尾部绝对误差界必须是有限非负数')
    entropy = math.log(2)-tail_integral+weights@energy
    output = dict(entropy_per_site=entropy.tolist(), weights=weights.tolist(),
                  statistical_covariance=None, statistical_stderr=None,
                  covariance_scope='energy covariance conditional on fixed tail_integral',
                  tail_energy_cross_covariance='not included; do not use for a jointly fitted tail',
                  tail_abs_bound=tail_abs_bound,
                  tail_interval=None if tail_abs_bound is None else
                  np.column_stack((entropy-tail_abs_bound, entropy+tail_abs_bound)).tolist(),
                  provisional=True, interpolation_uncertainty='unresolved',
                  scope='integration kernel only; not accepted lattice SSE or a plateau claim')
    if covariance is not None:
        covariance = np.asarray(covariance, dtype=float)
        if covariance.shape != weights.shape or not np.all(np.isfinite(covariance)):
            raise ValueError('协方差维数或数值错误')
        scale = float(np.max(np.abs(covariance)))
        tolerance = 64*np.finfo(float).eps*len(weights)*scale
        if np.max(np.abs(covariance-covariance.T)) > tolerance:
            raise ValueError('协方差不对称')
        covariance = (covariance+covariance.T)/2
        if np.min(np.linalg.eigvalsh(covariance)) < -tolerance:
            raise ValueError('协方差不是半正定矩阵')
        propagated = weights@covariance@weights.T
        output['statistical_covariance'] = propagated.tolist()
        output['statistical_stderr'] = np.sqrt(np.maximum(np.diag(propagated), 0)).tolist()
    return output


def self_check():
    from star_thermal_reference import reference

    beta_min, beta_max = 1/3, 100.
    hot = reference(1/beta_min)
    tail = math.log(2)+beta_min*hot['energy_per_site']-hot['entropy_per_site']
    grids = [np.array([1/3, 1., 1/.3, 10., 1/.03, 100.])]
    grids.extend(np.geomspace(beta_min, beta_max, n) for n in [33, 65, 129, 257])
    checks = []
    for beta in grids:
        exact = [reference(1/b) for b in beta]
        energy = [x['energy_per_site'] for x in exact]
        expected = np.array([x['entropy_per_site'] for x in exact])
        actual = reconstruct(beta, energy, tail)
        errors = np.array(actual['entropy_per_site'])-expected
        checks.append(dict(points=len(beta), max_abs_error=float(np.max(np.abs(errors))),
                           cold_entropy_error=float(errors[-1]),
                           uses_exact_tail=True, statistical_noise=False))
    refined = [x['max_abs_error'] for x in checks[1:]]
    passed = all(b < a for a, b in zip(refined, refined[1:])) and refined[-1] < 1e-4
    leading_tail = -9*beta_min**2/160
    paths = [Path(__file__), Path(__file__).with_name('star_thermal_reference.py'),
             Path(__file__).with_name('check_super_honeycomb_extensions.py')]
    return dict(status='passed' if passed else 'failed',
                scope='deterministic JB=0 quadrature check, not lattice SSE',
                criterion='33/65/129/257-point errors decrease; final max error < 1e-4',
                checks=checks, exact_tail=tail, leading_high_T_tail=leading_tail,
                entropy_bias_from_leading_tail=tail-leading_tail,
                tail_warning='leading high-T term has no supplied remainder bound',
                provenance={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-check', action='store_true', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = self_check()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    with (args.output_dir/'validation.json').open('x', encoding='utf-8') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result['status'] == 'passed' else 1)


if __name__ == '__main__':
    main()
