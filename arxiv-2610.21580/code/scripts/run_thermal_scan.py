"""运行有限规模热力学扫描，保留失败单元及原始分箱。"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def save(path, data):
    with path.open('x',encoding='utf-8') as handle:
        json.dump(data,handle,ensure_ascii=False,indent=2,allow_nan=False)
        handle.write('\n')


def mean_error(x):
    return statistics.fmean(x),statistics.stdev(x)/math.sqrt(len(x))


def lag1(x):
    m=statistics.fmean(x)
    d=sum((v-m)**2 for v in x)
    return sum((a-m)*(b-m) for a,b in zip(x,x[1:]))/d if d else 0.0


def heat_jackknife(n,n2,nsite):
    count=len(n)
    a,b=sum(n),sum(n2)
    full=(b/count-(a/count)**2-a/count)/nsite
    leave=[((b-y)/(count-1)-((a-x)/(count-1))**2-(a-x)/(count-1))/nsite for x,y in zip(n,n2)]
    center=statistics.fmean(leave)
    return count*full-(count-1)*center,math.sqrt((count-1)/count*sum((v-center)**2 for v in leave))


def analyze(path,nsite,beta,couplings,bond_counts):
    with path.open(encoding='utf-8') as handle:
        rows=list(csv.DictReader(handle))
    if len(rows)<16 or [int(r['bin']) for r in rows]!=list(range(1,len(rows)+1)):
        raise ValueError('分箱数量或编号错误')
    columns={k:[float(r[k]) for r in rows] for k in rows[0] if k!='bin'}
    if not all(math.isfinite(x) for c in columns.values() for x in c):
        raise ValueError('分箱含非有限值')
    n=columns['mean_expansion_order']
    n2=[c*nsite+x*x+x for c,x in zip(columns['specific_heat_per_site'],n)]
    if max(abs(a+b-c) for a,b,c in zip(columns['mean_order_type1'],columns['mean_order_type2'],n))>1e-8:
        raise ValueError('两类展开阶数之和不等于总阶数')
    quantities={}
    for k in ['energy_per_site','uniform_susceptibility','staggered_m2','ferrimagnetic_m2',
              'sublattice_plus_m2','sublattice_minus_m2','sublattice_cross']:
        value,error=mean_error(columns[k])
        quantities[k]=dict(value=value,error=error)
    heat,error=heat_jackknife(n,n2,nsite)
    quantities['specific_heat_per_site']=dict(value=heat,error=error,estimator='bias-corrected delete-one-bin jackknife on raw n,n2')
    value,error=mean_error([x/beta for x in columns['uniform_susceptibility']])
    quantities['T_chi']=dict(value=value,error=error)
    for t in [0,1]:
        corr=[0.25-x/(beta*bond_counts[t]*couplings[t]) for x in columns[f'mean_order_type{t+1}']]
        value,error=mean_error(corr)
        quantities[f'bond_correlation_type{t+1}']=dict(value=value,error=error)
    reconstructed=sum(quantities[f'bond_correlation_type{t+1}']['value']*bond_counts[t]*couplings[t]/nsite for t in [0,1])
    if abs(reconstructed-quantities['energy_per_site']['value'])>1e-10:
        raise ValueError('分键能量与总能量不一致')
    diagnostics={}
    for name in ['energy_per_site','uniform_susceptibility','mean_expansion_order']:
        x=columns[name]
        _,error=mean_error(x)
        blocked=[statistics.fmean(x[i:i+2]) for i in range(0,len(x)-1,2)]
        _,blocked_error=mean_error(blocked)
        half=len(x)//2
        a,ea=mean_error(x[:half]);b,eb=mean_error(x[half:])
        denominator=math.hypot(ea,eb)
        diagnostics[name]=dict(lag1=lag1(x),blocking2_ratio=blocked_error/error if error else 1.,
                                split_half_z=abs(a-b)/denominator if denominator else 0.)
    return quantities,diagnostics


def execute(cell,spec,root,executable):
    params=cell['params'];settings={**spec['settings'],**cell.get('settings',{})}
    graph=root/settings['graph']
    provenance=spec['provenance']
    for path,key in [(graph,'graph_sha256'),(executable,'executable_sha256'),
                     (root/'sse_graph_thermal_v3.f90','source_sha256'),(Path(__file__),'runner_sha256')]:
        if sha(path)!=provenance[key]: raise ValueError(f'来源哈希错误：{key}')
    cell_dir=root/spec['run_dir']/'cells'/cell['cell_id']
    cell_dir.mkdir(parents=True,exist_ok=False)
    size=int(params['L']);temperature=float(params['T']);couplings=[1.,float(params['ratio'])]
    if size%2 or size<2 or temperature<=0 or couplings[1]<=0: raise ValueError('参数非法')
    beta=1/temperature
    seed=int.from_bytes(hashlib.sha256((spec['run_id']+':'+cell['cell_id']).encode()).digest()[:4],'big')%2000000000+1
    text=f'{graph.resolve()}\n{size} {size}\n{beta:.17g} 2\n1 {couplings[1]:.17g}\n{settings["nbins"]} {settings["equil_sweeps"]} {settings["bin_sweeps"]} {seed}\n'
    with (cell_dir/'sse_graph.in').open('x',encoding='utf-8') as handle: handle.write(text)
    env={**os.environ,**{k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']}}
    start=time.monotonic();rc=125;stdout='';stderr=''
    try:
        if hasattr(os,'getloadavg'):
            if os.getloadavg()[1]/len(os.sched_getaffinity(0))>=0.80: raise RuntimeError('CPU 繁忙，停止该单元')
        run=subprocess.run([str(executable),str((cell_dir/'sse_graph.in').resolve())],cwd=cell_dir,
                           env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=settings['timeout_seconds'])
        rc,stdout,stderr=run.returncode,run.stdout,run.stderr
    except subprocess.TimeoutExpired:
        rc=124;stderr='单元超时，子进程已终止'
    except Exception as exc:
        stderr=str(exc)
    for name,content in [('stdout.txt',stdout),('stderr.txt',stderr)]:
        with (cell_dir/name).open('x',encoding='utf-8') as handle: handle.write(content)
    manifest=dict(schema_version=1,run_id=spec['run_id'],cell_id=cell['cell_id'],params=params,
                  settings=settings,provenance=provenance,status=dict(success=False),
                  execution=dict(returncode=rc,seed=seed,elapsed_seconds=time.monotonic()-start,beta=beta,couplings=couplings))
    if rc==0:
        try:
            observations,diagnostics=analyze(cell_dir/'sse_graph_bins.csv',settings['basis_sites']*size**2,beta,couplings,
                                            [x*size**2 for x in settings['bonds_per_type_cell']])
            manifest.update(observables=observations,diagnostics=diagnostics)
            manifest['status']['success']=True
        except Exception as exc:
            manifest['validation_error']=str(exc)
    manifest['artifact_sha256']={p.name:sha(p) for p in cell_dir.iterdir() if p.is_file()}
    save(cell_dir/'manifest.json',manifest)
    print(cell['cell_id'],manifest['status'],flush=True)
    return manifest['status']['success']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-spec',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    if not 1<=args.workers<=4: raise ValueError('并发上限为 4')
    spec=json.loads(args.run_spec.read_text(encoding='utf-8'))
    if len({c['cell_id'] for c in spec['cells']})!=len(spec['cells']): raise ValueError('重复单元')
    for cell in spec['cells']:
        if (args.root/spec['run_dir']/'cells'/cell['cell_id']).exists(): raise FileExistsError('已有单元，拒绝自动重投')
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(execute,c,spec,args.root.resolve(),args.executable.resolve()) for c in spec['cells']]
        successes=[f.result() for f in futures]
    save(args.root/spec['run_dir']/'batch_status.json',dict(run_id=spec['run_id'],completed=sum(successes),total=len(successes)))
    raise SystemExit(0 if all(successes) else 1)
