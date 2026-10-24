"""核验远端 CLS 谱矩的来源、参数和恒等式。"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def validate(stage, root):
    expected_names={'result.json','stdout.txt','stderr.txt','environment.txt','source.sha256','exit_code.txt'}
    entries={}
    for line in (stage/'artifacts.sha256').read_text().splitlines():
        checksum,name=line.split(maxsplit=1)
        require(name in expected_names and name not in entries,'工件名单错误')
        require(digest(stage/name)==checksum.lower(),f'工件哈希错误：{name}')
        entries[name]=checksum
    require(set(entries)==expected_names,'工件缺失')
    require((stage/'exit_code.txt').read_text().strip()=='0','运行返回码非零')
    require(not (stage/'stderr.txt').read_text().strip(),'运行含标准错误输出')
    sources={'ed_cls_moments.py':root/'scripts/ed_cls_moments.py',
             'ed_graph_heisenberg.py':root/'scripts/ed_graph_heisenberg.py',
             'super_honeycomb_poscar.dat':root/'graphs/super_honeycomb_poscar.dat'}
    source_entries={}
    for line in (stage/'source.sha256').read_text().splitlines():
        checksum,name=line.split(maxsplit=1)
        require(name in sources and name not in source_entries,'来源名单错误')
        require(digest(sources[name])==checksum.lower(),'远端来源哈希错误')
        source_entries[name]=checksum.lower()
    require(set(source_entries)==set(sources),'来源缺失')
    data=json.loads((stage/'result.json').read_text(encoding='utf-8'))
    def finite(value):
        if isinstance(value,float):
            require(math.isfinite(value),'含非有限数值')
        elif isinstance(value,dict):
            for item in value.values(): finite(item)
        elif isinstance(value,list):
            for item in value: finite(item)
    finite(data)
    for key,value in dict(schema_version=1,status='completed',lx=2,ly=2,nsite=20,nbond=24,
                          ground_spin=2,ground_magnetization=2).items():
        require(data[key]==value,f'参数错误：{key}')
    require(data['solver']==dict(method='eigsh',tol=1e-12,maxiter=1200,seed=20260917),'求解器参数错误')
    require(data['hamiltonian']=='H=sum_<ij> S_i dot S_j; S=1/2; J=1; h=0; PBC','模型错误')
    for key,name in [('source_sha256','ed_cls_moments.py'),('base_source_sha256','ed_graph_heisenberg.py'),
                     ('graph_sha256','super_honeycomb_poscar.dat')]:
        require(data[key].lower()==source_entries[name],'结果来源哈希错误')
    reference=json.loads((root/'results/ed_super_honeycomb_v1/result.json').read_text(encoding='utf-8'))
    require(reference['graph_sha256'].lower()==data['graph_sha256'].lower(),'基准图不一致')
    energy_difference=abs(data['ground_energy']-reference['ground_state_energy'])
    require(energy_difference<1e-9 and 0<=data['ground_residual']<1e-9,'基态能量或残差错误')
    weights=data['cls_weights']
    require(len(weights)==20 and sum(w!=0 for w in weights)==6,'不是六点支撑')
    require(abs(sum(w*w for w in weights)-1)<1e-12 and abs(sum(weights))<1e-12,'CLS 归一化错误')
    require(data['adjacency_cls_residual']<1e-12,'CLS 湮灭残差错误')
    controls=data['controls']
    for channel in ['minus','z']:
        require(abs(controls[channel]['raw_weight']-4)<1e-9,'总自旋权重错误')
        require(abs(controls[channel]['mean_energy'])<1e-9 and controls[channel]['energy_std']<1e-9,'总自旋零能检查失败')
    require(controls['plus']['raw_weight']<1e-18,'最高磁化态检查失败')
    probes={'local_cls','cls_q_1pi_0pi','cls_q_0pi_1pi','cls_q_1pi_1pi'}
    require(set(data['channels'])=={'minus','z','plus'},'通道缺失')
    for channel,items in data['channels'].items():
        require(set(items)==probes,'动量探针缺失')
        for name,item in items.items():
            require(not item['null_inelastic_state'],'探针意外为零')
            require(item['inelastic_weight']>0 and item['elastic_weight']>=0,'谱权重错误')
            require(abs(item['raw_weight']-item['elastic_weight']-item['inelastic_weight'])<1e-10,'权重分解错误')
            require(item['mean_energy']>=0 and item['energy_variance']>=0,'谱矩符号错误')
            require(abs(item['energy_std']**2-item['energy_variance'])<1e-10,'标准差错误')
            require(abs(item['second_moment']-item['mean_energy']**2-item['energy_variance'])<1e-10,'谱矩恒等式错误')
    environment=json.loads((stage/'environment.txt').read_text())
    require(environment['nproc']>0 and environment['load'][1]/environment['nproc']<0.80,'负载门槛失败')
    return dict(status='passed',nsite=20,probe_count=12,energy_reference_difference=energy_difference,
                result_sha256=digest(stage/'result.json'),validator_sha256=digest(Path(__file__)),
                evidence='File/source hashes, fixed parameters, reference energy, six-site CLS, SU(2) controls and moment identities checked; no independent N20 full-spectrum audit.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',type=Path)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    require(not args.report.exists(),'拒绝覆盖核验报告')
    result=validate(args.stage,args.root)
    with args.report.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,ensure_ascii=False,indent=2)
        handle.write('\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
