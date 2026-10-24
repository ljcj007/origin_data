"""将两类独立二聚体的 SSE 与精确热力学逐项比较。"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from run_thermal_scan import analyze,save,sha


ROOT = Path(__file__).resolve().parents[1]
SETTINGS = dict(graph='graphs/two_coupling_dimers.dat', basis_sites=4,
                bonds_per_type_cell=[1, 1], nbins=64, equil_sweeps=5000,
                bin_sweeps=1000, timeout_seconds=180)
POINTS = set(itertools.product([2], [.1, 1.], [.05, .3, 2.], [1, 2]))
REQUIRED = {'sse_graph.in', 'sse_graph_bins.csv', 'sse_graph_summary.txt',
            'stdout.txt', 'stderr.txt'}
BIN_COLUMNS = set(('bin energy_per_site specific_heat_per_site staggered_m2 staggered_m4 '
                   'binder_z_staggered ferrimagnetic_m2 uniform_susceptibility '
                   'sublattice_plus_m2 sublattice_minus_m2 sublattice_cross '
                   'spin_stiffness_cell_x spin_stiffness_cell_y spin_stiffness_cell_average '
                   'spin_stiffness_times_L winding_squared mean_expansion_order '
                   'mean_order_type1 mean_order_type2').split())


def close(actual, expected, label):
    if not isinstance(actual, (int, float)) or not math.isfinite(actual):
        raise ValueError(f'非有限数值：{label}')
    if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-10):
        raise ValueError(f'数值不一致：{label}')


def check_spec(spec, batch, run_dir):
    run_id = spec['run_id']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', run_id) or run_dir.name != run_id:
        raise ValueError('批次目录与编号不一致')
    if spec['run_dir'] != f'results/{run_id}' or spec['settings'] != SETTINGS:
        raise ValueError('批次路径或参数偏离已登记二聚体协议')
    cells = spec['cells']
    ids = [c['cell_id'] for c in cells]
    if len(cells) != 12 or set(ids) != {f'cell-{i:04d}' for i in range(1, 13)}:
        raise ValueError('必须恰有十二个不同单元')
    points = []
    for cell in cells:
        params = cell['params']
        if set(params) != {'L', 'ratio', 'T', 'replica'} or cell.get('settings'):
            raise ValueError('单元参数或覆盖项错误')
        points.append(tuple(params[k] for k in ['L', 'ratio', 'T', 'replica']))
    if set(points) != POINTS or len(set(points)) != 12:
        raise ValueError('温度、交换比或独立副本不完整')
    if batch != dict(run_id=run_id, completed=12, total=12):
        raise ValueError('批次未全部成功完成')
    if {p.name for p in (run_dir/'cells').iterdir()} != set(ids):
        raise ValueError('实际单元目录与计划不一致')


def check_input(path, params, seed):
    rows = path.read_text(encoding='utf-8').splitlines()
    if len(rows) != 5 or not rows[0].replace('\\', '/').endswith('/'+SETTINGS['graph']):
        raise ValueError('输入图路径或行数错误')
    if list(map(int, rows[1].split())) != [params['L']]*2:
        raise ValueError('输入尺寸错误')
    beta, types = rows[2].split()
    close(float(beta), 1/params['T'], '输入 beta')
    if int(types) != 2 or list(map(float, rows[3].split())) != [1., params['ratio']]:
        raise ValueError('输入交换常数错误')
    if list(map(int, rows[4].split())) != [64, 5000, 1000, seed]:
        raise ValueError('输入分箱、步数或种子错误')


def check_bins(path, params):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        if (set(reader.fieldnames or []) != BIN_COLUMNS
                or len(reader.fieldnames) != len(BIN_COLUMNS)):
            raise ValueError('分箱必须包含箱号及全部十八个观测量，且不能重列')
        rows = list(reader)
    if len(rows) != SETTINGS['nbins']:
        raise ValueError('原始分箱不是已登记的 64 箱')
    beta, nsite = 1/params['T'], 4*params['L']**2
    offset = (1+params['ratio'])/16
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError('分箱列数不一致')
        values = {k: float(v) for k, v in row.items() if k != 'bin'}
        if not all(math.isfinite(v) for v in values.values()):
            raise ValueError('分箱包含非有限数')
        close(values['energy_per_site'], offset-values['mean_expansion_order']/(beta*nsite), '逐箱能量')
        close(values['uniform_susceptibility'], beta*nsite*values['ferrimagnetic_m2']/3, '逐箱磁化率')
        plus, minus, cross = [values[k] for k in ['sublattice_plus_m2', 'sublattice_minus_m2', 'sublattice_cross']]
        close(values['staggered_m2'], (plus+minus-2*cross)/4, '逐箱交错关联')
        close(values['ferrimagnetic_m2'], (plus+minus+2*cross)/4, '逐箱均匀关联')
        for name in ['spin_stiffness_cell_x', 'spin_stiffness_cell_y', 'spin_stiffness_cell_average',
                     'spin_stiffness_times_L', 'winding_squared']:
            close(values[name], 0., name)


def check_cell(cell, spec, run_dir):
    folder = run_dir/'cells'/cell['cell_id']
    if folder.resolve().parent != (run_dir/'cells').resolve():
        raise ValueError('单元目录越界')
    m = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    if (m['run_id'] != spec['run_id'] or m['cell_id'] != cell['cell_id']
            or m['status']['success'] is not True or m['execution']['returncode'] != 0
            or m['params'] != cell['params'] or m['settings'] != SETTINGS
            or m['provenance'] != spec['provenance']):
        raise ValueError('运行清单或返回码不一致')
    artifacts = m['artifact_sha256']
    actual = {p.name for p in folder.iterdir() if p.is_file() and p.name != 'manifest.json'}
    if not REQUIRED.issubset(artifacts) or set(artifacts) != actual:
        raise ValueError('工件缺失或未登记')
    for name, digest in artifacts.items():
        path = folder/name
        if Path(name).name != name or path.resolve().parent != folder.resolve() or sha(path) != digest:
            raise ValueError('工件路径或哈希错误')
    seed = int.from_bytes(hashlib.sha256((spec['run_id']+':'+cell['cell_id']).encode()).digest()[:4], 'big')%2000000000+1
    if m['execution']['seed'] != seed or m['execution']['couplings'] != [1., cell['params']['ratio']]:
        raise ValueError('种子或交换记录错误')
    close(m['execution']['beta'], 1/cell['params']['T'], '执行 beta')
    check_input(folder/'sse_graph.in', cell['params'], seed)
    check_bins(folder/'sse_graph_bins.csv', cell['params'])
    l = cell['params']['L']
    observations, diagnostics = analyze(folder/'sse_graph_bins.csv', 4*l*l,
                                       1/cell['params']['T'], [1., cell['params']['ratio']], [l*l]*2)
    if set(observations) != set(m['observables']) or set(diagnostics) != set(m['diagnostics']):
        raise ValueError('汇总观测量或诊断项缺失')
    for label, recalculated, recorded in [('观测量', observations, m['observables']),
                                         ('诊断', diagnostics, m['diagnostics'])]:
        for name, fields in recalculated.items():
            if set(fields) != set(recorded[name]):
                raise ValueError('汇总字段缺失')
            for key, value in fields.items():
                if isinstance(value, str):
                    if value != recorded[name][key]:
                        raise ValueError('估计量描述错误')
                else:
                    close(recorded[name][key], value, f'{label}:{name}:{key}')
    return observations


def exact(ratio,temperature):
    correlations=[];heats=[];tchis=[]
    for coupling in [1.,ratio]:
        x=coupling/temperature;q=math.exp(-x);denominator=1+3*q
        correlations.append((-0.75+0.75*q)/denominator)
        heats.append(3*x*x*q/denominator**2)
        tchis.append(2*q/denominator)
    return dict(energy_per_site=(correlations[0]+ratio*correlations[1])/4,
                specific_heat_per_site=sum(heats)/4,T_chi=sum(tchis)/4,
                bond_correlation_type1=correlations[0],bond_correlation_type2=correlations[1])


def validate(run_dir, executable=None, root=ROOT):
    run_dir, root = Path(run_dir).resolve(), Path(root).resolve()
    spec=json.loads((run_dir/'run_spec.json').read_text(encoding='utf-8'))
    batch=json.loads((run_dir/'batch_status.json').read_text(encoding='utf-8'))
    check_spec(spec, batch, run_dir)
    if executable is None:
        raise ValueError('必须提供实际使用的可执行文件用于来源核对')
    graph = root/SETTINGS['graph']
    if graph.read_text().split() != '4 2 3 1 2 0 0 0 1 3 4 0 0 0 2 1 -1 1 -1 0 0 0'.split():
        raise ValueError('图不是已登记的两类独立二聚体')
    for path, key in [(graph, 'graph_sha256'), (Path(executable), 'executable_sha256'),
                      (root/'sse_graph_thermal_v3.f90', 'source_sha256'),
                      (root/'scripts/run_thermal_scan.py', 'runner_sha256')]:
        if sha(path) != spec['provenance'][key]:
            raise ValueError(f'实际来源哈希不一致：{key}')
    reports=[]
    for cell in spec['cells']:
        observations=check_cell(cell,spec,run_dir)
        oracle=exact(cell['params']['ratio'],cell['params']['T'])
        for name,expected in oracle.items():
            obs=observations[name];delta=obs['value']-expected
            passed=abs(delta)<=5*obs['error']+1e-4
            reports.append(dict(cell_id=cell['cell_id'],observable=name,exact=expected,**obs,
                                delta=delta,passed=passed))
    return dict(status='passed' if all(r['passed'] for r in reports) else 'failed',
                comparison_count=len(reports),checks=reports,validator_sha256=sha(Path(__file__)),
                run_id=spec['run_id'],run_spec_sha256=sha(run_dir/'run_spec.json'),
                integrity='12 cells, 64 bins/cell, source and artifact hashes, raw-bin recomputation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run_dir',type=Path)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--root',type=Path,default=ROOT)
    args=p.parse_args();result=validate(args.run_dir,args.executable,args.root)
    save(args.run_dir/'exact_validation.json',result)
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['status']=='passed' else 1)
