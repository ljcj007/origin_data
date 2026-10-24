"""从已下载原始分箱核对试扫描，保留全部统计提示。"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path, PurePosixPath
from types import FunctionType

from run_thermal_scan import analyze, heat_jackknife, save, sha
from validate_thermal_dimers import BIN_COLUMNS, REQUIRED, close, validate


ROOT = Path(__file__).resolve().parents[1]


def left_to_right_sum(values):
    """复现服务器 Python 3.10 的浮点求和顺序。"""
    total = 0
    for value in values:
        total += value
    return total


SERVER_HEAT = FunctionType(heat_jackknife.__code__,
                          {**heat_jackknife.__globals__, 'sum': left_to_right_sum})


def audit(root):
    root = Path(root).resolve()
    executable = root/'build/thermal_v3/sse_thermal'
    benchmark = validate(root/'results/thermal_dimers_remote_v1', executable, root)
    if benchmark['status'] != 'passed':
        raise ValueError('基准未通过')
    folder = root/'results/thermal_pilot_remote_v1'
    spec = json.loads((folder/'run_spec.json').read_text())
    expected = json.loads((ROOT/'results/thermal_pilot_remote_v1/run_spec.json').read_text())
    if spec != expected or len(spec['cells']) != 72:
        raise ValueError('运行清单与本地登记不一致')
    points = {tuple(c['params'][k] for k in ['L', 'ratio', 'T', 'replica']) for c in spec['cells']}
    if points != set(itertools.product([4, 8], [.02, .1, 1.], [.01, .03, .1, .3, 1., 3.], [1, 2])):
        raise ValueError('试扫描网格错误')
    settings = dict(graph='graphs/super_honeycomb_JA_JB.dat', basis_sites=5,
                    bonds_per_type_cell=[3, 3], nbins=32, equil_sweeps=5000,
                    bin_sweeps=500, timeout_seconds=300)
    if spec['settings'] != settings:
        raise ValueError('采样参数错误')
    sources = ['sse_graph_thermal_v3.f90', 'scripts/run_thermal_scan.py',
               'scripts/compile_thermal.py', 'scripts/validate_thermal_dimers.py',
               'graphs/super_honeycomb_JA_JB.dat', 'graphs/two_coupling_dimers.dat']
    for source in sources:
        if sha(root/source) != sha(ROOT/source):
            raise ValueError('下载来源与本地不一致：'+source)
    for path, key in [(root/sources[0], 'source_sha256'), (root/sources[1], 'runner_sha256'),
                      (root/settings['graph'], 'graph_sha256'), (executable, 'executable_sha256')]:
        if sha(path) != spec['provenance'][key]:
            raise ValueError('来源哈希错误：'+key)
    ids = [c['cell_id'] for c in spec['cells']]
    if (len(set(ids)) != 72 or set(ids) != {f'cell-{i:04d}' for i in range(1, 73)}
            or {p.name for p in (folder/'cells').iterdir()} != set(ids)):
        raise ValueError('实际单元与登记不一致')
    records, flags, successful, roundoff = [], [], 0, []
    for cell in spec['cells']:
        path = folder/'cells'/cell['cell_id']
        m = json.loads((path/'manifest.json').read_text())
        params = cell['params']
        for key, value in dict(run_id=spec['run_id'], cell_id=cell['cell_id'], params=params,
                               settings=settings, provenance=spec['provenance']).items():
            if m[key] != value:
                raise ValueError('单元记录不一致：'+cell['cell_id']+':'+key)
        files = {p.name for p in path.iterdir() if p.is_file() and p.name != 'manifest.json'}
        if files != set(m['artifact_sha256']):
            raise ValueError('工件覆盖不完整')
        for name, digest in m['artifact_sha256'].items():
            if Path(name).name != name or (path/name).resolve().parent != path.resolve() or sha(path/name) != digest:
                raise ValueError('工件路径或哈希错误')
        rc = m['execution']['returncode']
        record = dict(cell_id=cell['cell_id'], **params, success=m['status']['success'], returncode=rc)
        if m['status']['success'] is not True or rc != 0:
            record['failure'] = m.get('validation_error', (path/'stderr.txt').read_text())
            records.append(record)
            continue
        if not REQUIRED.issubset(files):
            raise ValueError('成功单元缺少工件')
        beta, size = 1/params['T'], params['L']
        seed = int.from_bytes(hashlib.sha256((spec['run_id']+':'+cell['cell_id']).encode()).digest()[:4], 'big')%2000000000+1
        if m['execution']['seed'] != seed or m['execution']['couplings'] != [1., params['ratio']]:
            raise ValueError('随机种子或交换记录错误')
        close(m['execution']['beta'], beta, 'beta')
        lines = (path/'sse_graph.in').read_text().splitlines()
        if (len(lines) != 5 or lines[0] != str(PurePosixPath(spec['provenance']['remote_root'])/settings['graph'])
                or list(map(int, lines[1].split())) != [size, size]
                or list(map(float, lines[2].split())) != [beta, 2.]
                or list(map(float, lines[3].split())) != [1., params['ratio']]
                or list(map(int, lines[4].split())) != [32, 5000, 500, seed]):
            raise ValueError('实际输入与计划不一致')
        with (path/'sse_graph_bins.csv').open(newline='') as handle:
            reader = csv.DictReader(handle)
            if set(reader.fieldnames or []) != BIN_COLUMNS or len(reader.fieldnames) != len(BIN_COLUMNS):
                raise ValueError('分箱列不完整')
            bins = list(reader)
        if len(bins) != 32:
            raise ValueError('成功单元不是32箱')
        for row in bins:
            if None in row or any(v is None for v in row.values()):
                raise ValueError('分箱列数错误')
            r = {k: float(v) for k, v in row.items()}
            nsite = 5*size*size
            close(r['energy_per_site'], 3*(1+params['ratio'])/20-r['mean_expansion_order']/(beta*nsite), '逐箱能量')
            close(r['uniform_susceptibility'], beta*nsite*r['ferrimagnetic_m2']/3, '逐箱磁化率')
            plus, minus, cross = (r[k] for k in ['sublattice_plus_m2', 'sublattice_minus_m2', 'sublattice_cross'])
            close(r['staggered_m2'], .16*plus+.36*minus-.48*cross, '逐箱交错关联')
            close(r['ferrimagnetic_m2'], .16*plus+.36*minus+.48*cross, '逐箱均匀关联')
            close(r['spin_stiffness_cell_average'], (r['spin_stiffness_cell_x']+r['spin_stiffness_cell_y'])/2, '刚度方向平均')
            close(r['spin_stiffness_times_L'], r['spin_stiffness_cell_average']*size, '刚度尺度')
            close(r['spin_stiffness_cell_average'], .75*r['winding_squared']/beta, '刚度绕数')
        observations, diagnostics = analyze(path/'sse_graph_bins.csv', 5*size*size, beta,
                                             [1., params['ratio']], [3*size*size]*2)
        n = [float(r['mean_expansion_order']) for r in bins]
        n2 = [float(r['specific_heat_per_site'])*nsite+x*x+x for r, x in zip(bins, n)]
        server_value, server_error = SERVER_HEAT(n, n2, nsite)
        replay = {**observations, 'specific_heat_per_site':
                  {**observations['specific_heat_per_site'], 'value': server_value, 'error': server_error}}
        roundoff.append(abs(observations['specific_heat_per_site']['value']-server_value))
        for label, actual in [('observables', replay), ('diagnostics', diagnostics)]:
            if actual.keys() != m[label].keys():
                raise ValueError('汇总字段不完整')
            for name, fields in actual.items():
                if fields.keys() != m[label][name].keys():
                    raise ValueError('汇总子字段不完整')
                for key, value in fields.items():
                    if isinstance(value, str):
                        if value != m[label][name][key]:
                            raise ValueError('估计量描述不一致')
                    else:
                        close(m[label][name][key], value, name+':'+key)
        for name, d in diagnostics.items():
            for key, limit in [('lag1', .3), ('blocking2_ratio', 1.5), ('split_half_z', 3.)]:
                if abs(d[key]) > limit:
                    flags.append(dict(cell_id=cell['cell_id'], observable=name, diagnostic=key, value=d[key]))
        record.update(observables=observations, diagnostics=diagnostics,
                      elapsed_seconds=m['execution']['elapsed_seconds'], bins=32)
        records.append(record)
        successful += 1
    batch = json.loads((folder/'batch_status.json').read_text())
    if batch != dict(run_id=spec['run_id'], completed=successful, total=72):
        raise ValueError('批次成功数与单元不一致')
    pairs = []
    grouped = {}
    for r in records:
        grouped.setdefault((r['L'], r['ratio'], r['T']), []).append(r)
    for point, replicas in grouped.items():
        if len(replicas) != 2 or {r['replica'] for r in replicas} != {1, 2}:
            raise ValueError('独立副本覆盖错误')
        if not all(r['success'] for r in replicas):
            continue
        a, b = replicas
        for name in a['observables']:
            x, y = a['observables'][name], b['observables'][name]
            error = math.hypot(x['error'], y['error'])
            delta = x['value']-y['value']
            z = abs(delta)/error if error else (0. if delta == 0 else None)
            pairs.append(dict(L=point[0], ratio=point[1], T=point[2], observable=name,
                              delta=delta, combined_error=error, z=z, flagged=z is None or z > 3))
    return dict(run_id=spec['run_id'], integrity='passed', completed=successful, total=72,
                raw_bins=successful*32, benchmark_checks=benchmark['comparison_count'],
                records=records, statistical_flags=flags, replica_comparisons=pairs,
                provisional=True, interpretation='coarse pilot only; no entropy plateau or phase-transition claim',
                floating_point_replay=dict(heat_sum_order='Python 3.10 left-to-right',
                                           max_center_difference=max(roundoff, default=0.),
                                           comparison_tolerance_unchanged=True),
                provenance=dict(archive_root=str(root), audit_sha256=sha(Path(__file__)),
                                run_spec_sha256=sha(folder/'run_spec.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.root)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    save(args.output_dir/'audit.json', result)
    rows = []
    for record in result['records']:
        row = {k: record[k] for k in ['cell_id', 'L', 'ratio', 'T', 'replica', 'success', 'returncode']}
        for name, item in record.get('observables', {}).items():
            row[name] = item['value']
            row[name+'_stderr'] = item['error']
        rows.append(row)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (args.output_dir/'observables.csv').open('x', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({k: result[k] for k in ['integrity', 'completed', 'total', 'raw_bins', 'provisional']}))
    print('统计提示', len(result['statistical_flags']), '副本提示', sum(x['flagged'] for x in result['replica_comparisons']))
