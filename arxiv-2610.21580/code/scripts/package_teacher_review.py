"""按白名单生成教师审阅包，逐文件校验且不覆盖已有压缩包。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_FILES = (
    'POSCAR (2).vasp', 'sse_graph.f90', 'sse_graph_stiffness_v2.f90',
    'sse_graph_thermal_v3.f90', 'ising_graph.f90',
    'paper/manuscript_v2.md', 'paper/manuscript_v2.tex', 'paper/manuscript_v2.pdf',
    'paper/manuscript_v2.bbl', 'paper/references_v2.bib',
    'paper/manuscript_v2_thermal_backfill_notes.md',
    'figures/figure1_super_honeycomb.pdf', 'figures/figure2_ising.pdf',
    'figures/figure3_sse.pdf', 'output/pdf/teacher_brief_20260906.pdf',
    'docs/teacher_review_20260906.md', 'docs/next_steps_20260906.md',
    'docs/thermal_run_results_20260906.md', 'docs/entropy_integration_audit_20260906.md',
    'docs/closest_spin_literature_check_20260906.md',
    'docs/super_honeycomb_model.md', 'docs/stiffness_normalization.md',
    'docs/sse_stationarity_audit.md', 'docs/fortran_sse_tutorial.md',
    'docs/lattice_graph_format.md', 'docs/lswt_super_honeycomb.md',
    'docs/lswt_low_temperature.md', 'docs/flatband_cls_topology.md',
    'docs/high_field_extension.md', 'docs/ising_estimator_boundary.md',
    'docs/super_honeycomb_N20_findings_and_research_scope_20260906.md',
    'docs/extension_completion_checklist.json',
    'results/entropy_integration_audit_20260906.json',
    'results/literature_graph_comparison_20260906.json',
    'results/thermal_remote_execution_20260906.json',
    'results/manuscript_v2_thermal_audit_20260906.json',
    'results/teacher_backfill_verification_20260906.json',
)
SCRIPT_NAMES = (
    'analyze_ising_pilot.py', 'analyze_ising_scan.py', 'analyze_ising_universality.py',
    'analyze_sse_beta.py', 'analyze_sse_fss.py', 'audit_cls_moments_small.py',
    'audit_ising_wolff.py', 'audit_manuscript.py', 'audit_sse_stationarity.py',
    'audit_thermal_pilot.py', 'band_graph.py', 'build_manuscript_v2.ps1',
    'build_teacher_brief.py',
    'check_super_honeycomb_extensions.py', 'checksum_ledger.py',
    'compare_literature_graphs.py', 'compile_thermal.py', 'ed_cls_moments.py',
    'ed_graph_heisenberg.py', 'ed_square_octagon.py', 'ed_super_honeycomb_xdiag.jl',
    'figure1_super_honeycomb.py', 'figure2_ising.py', 'figure3_sse.py',
    'flatband_cls_topology.py', 'ising_l2_exact.py', 'lattice_from_poscar.py',
    'lswt_low_temperature.py', 'lswt_super_honeycomb.py', 'one_magnon_high_field.py',
    'package_teacher_review.py', 'parameter_scan_cli.py', 'peripheral_benchmarks.py',
    'run_ising_cell.py', 'run_sse_cell.py', 'run_sse_cell_v3.py', 'run_thermal_scan.py',
    'star_thermal_reference.py', 'summarize_ising_universality.py',
    'summarize_sse_fss_fits.py', 'thermal_entropy_integrator.py',
    'validate_cls_moments.py', 'validate_ed_super_honeycomb.py',
    'validate_ising_l2_exact.py', 'validate_sse_dimer.py', 'validate_thermal_dimers.py',
)
TEST_NAMES = (
    'analyze_ising_universality', 'analyze_sse_fss', 'audit_ising_wolff',
    'audit_manuscript', 'audit_sse_stationarity', 'audit_thermal_pilot',
    'band_graph', 'checksum_ledger', 'compare_literature_graphs', 'ed_cls_moments',
    'figure3_sse', 'flatband_cls_topology', 'ising_l2_exact', 'lattice_from_poscar',
    'lswt_low_temperature', 'lswt_super_honeycomb', 'manuscript_v2_gauge',
    'manuscript_v2_literature', 'manuscript_v2_thermal', 'one_magnon_high_field',
    'package_teacher_review', 'peripheral_benchmarks', 'star_thermal_reference', 'thermal_entropy_integrator',
    'thermal_scan', 'validate_cls_moments', 'validate_ed_super_honeycomb',
    'validate_sse_dimer', 'validate_thermal_dimers',
)
DIRECTORIES = (
    'graphs', 'tests/fixtures', 'analysis/thermal_v1',
    'results/thermal_remote_complete_v1', 'results/thermal_pilot_audit_v1',
    'results/thermal_dimers_remote_v1', 'results/thermal_pilot_remote_v1',
    'results/sse_fss_v1', 'results/sse_fss_stationarity_v1', 'results/sse_beta_pilot_v2',
    'results/ed_super_honeycomb_v1', 'results/cls_moments_N10_20260905_v1',
    'results/cls_moments_N20_20260906_v1', 'results/flatband_cls_topology_v1',
    'results/lswt_super_honeycomb_v1', 'results/lswt_low_temperature_v1',
    'results/one_magnon_high_field_v1', 'results/star_thermal_reference_v1',
    'results/entropy_quadrature_reference_v2',
    'results/super_honeycomb_literature_extension_v1', 'results/literature_search_v1',
    'results/ising_tc_wolff_v1', 'results/ising_tc_wolff_audit_v1',
    'results/ising_tc_wolff_scaling_v1', 'results/ising_l2_exact_validation_v1',
    'results/peripheral_benchmarks_v1', 'results/sse_dimer_validation_v1',
    'results/poscar2',
)
SEALED_BUILD_PREFIXES = (
    'results/thermal_remote_complete_v1/',
    'results/peripheral_benchmarks_v1/', 'results/sse_dimer_validation_v1/',
)
FORBIDDEN_PARTS = {'.ssh', '.remote', '.git', '.tools', 'remote-bin', 'node_modules'}
FORBIDDEN_SUFFIXES = {'.pem', '.key', '.pfx', '.p12', '.ppk', '.kdbx', '.env'}
FORBIDDEN_NAMES = {'.env', 'id_rsa', 'id_ed25519', 'id_ecdsa', 'known_hosts', 'authorized_keys'}
SECRET_PATTERNS = (
    re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
    re.compile(rb'gh[pousr]_[A-Za-z0-9]{30,}'),
    re.compile(rb'github_pat_[A-Za-z0-9_]{40,}'),
    re.compile(rb'AKIA[0-9A-Z]{16}'),
)
MAX_FILE_BYTES = 8 * 1024 * 1024
MAX_ARCHIVE_BYTES = 32 * 1024 * 1024


def digest(data):
    return hashlib.sha256(data).hexdigest().upper()


def safe_name(name):
    path = PurePosixPath(name)
    if (not name or '\\' in name or ':' in name or path.is_absolute()
            or any(p in ('', '.', '..') for p in name.split('/'))
            or str(path) != name or any(p.lower() in FORBIDDEN_PARTS for p in path.parts)
            or path.suffix.lower() in FORBIDDEN_SUFFIXES
            or path.name.lower() in FORBIDDEN_NAMES or path.name.lower().startswith('.env.')):
        raise ValueError('不安全的包内路径：' + name)
    sealed = name.startswith(SEALED_BUILD_PREFIXES)
    cache = '__pycache__' in path.parts or path.suffix.lower() == '.pyc'
    if cache and not name.startswith('results/thermal_remote_complete_v1/'):
        raise ValueError('未登记的缓存：' + name)
    if path.suffix.lower() in {'.mod', '.exe'} and not sealed:
        raise ValueError('未登记的缓存或二进制：' + name)


def safe_content(name, data):
    safe_name(name)
    if len(data) > MAX_FILE_BYTES:
        raise ValueError('文件超过大小上限：' + name)
    if any(pattern.search(data) for pattern in SECRET_PATTERNS):
        raise ValueError('疑似凭据，停止打包：' + name)


def collect(root):
    names = set(REQUIRED_FILES)
    names.update('scripts/' + name for name in SCRIPT_NAMES)
    names.update('tests/test_' + name + '.py' for name in TEST_NAMES)
    missing = []
    for directory in DIRECTORIES:
        folder = root / directory
        if not folder.is_dir():
            missing.append(directory + '/')
            continue
        for path in folder.rglob('*'):
            if path.is_symlink():
                raise ValueError('不打包符号链接：' + str(path))
            if path.is_file():
                names.add(path.relative_to(root).as_posix())
    entries = []
    for name in sorted(names):
        safe_name(name)
        path = root / name
        if not path.is_file():
            missing.append(name)
            continue
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError('路径越界：' + name)
        data = path.read_bytes()
        safe_content(name, data)
        entries.append(dict(path=name, size=len(data), sha256=digest(data)))
    return entries, sorted(missing)


def scope_text(entries):
    return '''# 教师审阅包范围

这是按白名单整理的研究审阅与数据复核包，不是整个工作区备份，也不是正式投稿包。

先读 `docs/teacher_review_20260906.md` 或 `output/pdf/teacher_brief_20260906.pdf`；
论文为 `paper/manuscript_v2.pdf`，后续计划为 `docs/next_steps_20260906.md`。
两交换热扫描仍是试扫描，保留全部统计提示；不据此声称熵平台、双比热峰或有限温相变。

## 包含

- 当前论文 Markdown、REVTeX、参考文献、PDF，以及三幅 PDF 图。
- 核心物理程序、离线分析脚本、精选测试、晶格和 POSCAR。
- SSE 有限尺寸分箱、已采用的 beta 收敛扫描、Ising Wolff 数据及审计。
- ED、紧局域激发检验、星团、低温理论、文献图比较与熵求积验证结果。
- 12 单元热基准和 72 单元试扫描的原始分箱、清单、来源及工件账本。

## 明确不包含

- SSH 密钥、连接配置、远端控制或上传助手、系统环境、安装缓存。
- 失败或失效的旧试扫描、重复下载暂存区、整个历史工作区。
- 他人论文原 PDF、高分辨率 TIFF、与本文无关的临时文件。
- 封存 v1 全套论文；当前论文沿用的数值证据已单独收录。

为不破坏既有账本，三个封存验证目录中的编译工件原样保留；
thermal 封存目录还含两个当时生成的 Python 缓存。它们仅供哈希核验，不应直接执行。
历史运行记录中的本地和服务器路径是来源元数据，不是可用的连接配置。
本包不提供服务器登录能力，不会自动提交计算，也不自动安装依赖。

## 离线复核

以下命令从解压目录执行；使用已有 Python 3.10 或更新版本。
热扫描审计只用标准库，不执行封存的 Linux 二进制。

```text
python scripts/checksum_ledger.py results/thermal_remote_complete_v1 --verify
python scripts/audit_thermal_pilot.py --root results/thermal_remote_complete_v1 --output-dir review_thermal_audit
python -m unittest discover -s tests -p test_audit_thermal_pilot.py -v
python -m unittest discover -s tests -p test_validate_thermal_dimers.py -v
python scripts/package_teacher_review.py --verify ../teacher_review_20260906.zip
```

审计输出目录须尚不存在。预期为 72/72、2304 个原始分箱、provisional=true；
26 条单元统计提示与 2 条副本提示应保留，不作为删点依据。
其他分析可能需要已有 NumPy、SciPy、Matplotlib；编译轨迹测试需要已有 gfortran，
缺少编译器时该测试会显式跳过。论文构建需要已有 REVTeX、BibTeX 与 LaTeX，
`scripts/build_manuscript_v2.ps1` 的选项针对 MiKTeX；TeX Live 用户应按本机选项构建。
所附测试是离线子集，不代表完整工作区的服务器安全测试已打包。

`PACKAGE_MANIFEST.json` 为本次包内逐文件 SHA-256 账本（不包含它自身）；
历史账本是原始来源记录，里面的绝对路径不会因解压位置而改写。
''' + '\n本次原始文件数：' + str(len(entries)) + '。\n'


def verify_thermal_ledger(root, entries):
    prefix = 'results/thermal_remote_complete_v1/'
    ledger = json.loads((root / prefix / 'checksums.sha256.json').read_text(encoding='utf-8'))
    actual = [dict(path=e['path'][len(prefix):], size=e['size'], sha256=e['sha256'])
              for e in entries if e['path'].startswith(prefix)
              and e['path'] != prefix + 'checksums.sha256.json']
    if ledger['file_count'] != 522 or len(actual) != 522 or ledger['files'] != actual:
        raise ValueError('封存 thermal 的 522 文件账本不完整或不一致')


def verify(archive):
    with zipfile.ZipFile(archive) as handle:
        infos = handle.infolist()
        names = [i.filename for i in infos]
        if len(names) != len(set(names)) or len({n.casefold() for n in names}) != len(names):
            raise ValueError('重复或大小写冲突的 ZIP 路径')
        if sum(i.file_size for i in infos) > MAX_ARCHIVE_BYTES:
            raise ValueError('解压大小超限')
        for info in infos:
            safe_name(info.filename)
            if info.is_dir() or stat.S_ISLNK(info.external_attr >> 16) or info.flag_bits & 1:
                raise ValueError('不接受目录项、链接或加密项')
            if info.file_size > MAX_FILE_BYTES:
                raise ValueError('ZIP 单文件大小超限')
        manifest = json.loads(handle.read('PACKAGE_MANIFEST.json'))
        entries = manifest['files']
        if (manifest['file_count'] != len(entries)
                or len({e['path'] for e in entries}) != len(entries)
                or set(names) != {e['path'] for e in entries} | {'PACKAGE_MANIFEST.json'}):
            raise ValueError('ZIP 与清单文件集不一致')
        for entry in entries:
            data = handle.read(entry['path'])
            safe_content(entry['path'], data)
            if len(data) != entry['size'] or digest(data) != entry['sha256']:
                raise ValueError('ZIP 文件哈希不一致：' + entry['path'])
    return dict(status='verified', archive_files=len(names),
                uncompressed_bytes=sum(i.file_size for i in infos),
                archive_bytes=archive.stat().st_size, sha256=digest(archive.read_bytes()))


def package(root, output, entries):
    if output.exists():
        raise FileExistsError('不覆盖已有压缩包：' + str(output))
    payload = {}
    for entry in entries:
        data = (root / entry['path']).read_bytes()
        if len(data) != entry['size'] or digest(data) != entry['sha256']:
            raise ValueError('文件在盘点后变化：' + entry['path'])
        payload[entry['path']] = data
    payload['PACKAGE_SCOPE.md'] = scope_text(entries).encode('utf-8')
    records = [dict(path=name, size=len(data), sha256=digest(data))
               for name, data in sorted(payload.items())]
    manifest = dict(schema_version=1, scope='allowlisted teacher review; not a full backup',
                    file_count=len(records), files=records)
    payload['PACKAGE_MANIFEST.json'] = (json.dumps(manifest, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
    if sum(map(len, payload.values())) > MAX_ARCHIVE_BYTES:
        raise ValueError('包内容总量超限')
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.teacher-review-', suffix='.zip', dir=output.parent)
    os.close(fd)
    staged = Path(temporary)
    try:
        with zipfile.ZipFile(staged, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
            for name, data in sorted(payload.items()):
                safe_content(name, data)
                info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = (stat.S_IFREG | 0o644) << 16
                handle.writestr(info, data, compresslevel=9)
        verify(staged)
        # 硬链接发布，已有目标不会被覆盖。
        os.link(staged, output)
    finally:
        staged.unlink(missing_ok=True)
    return verify(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=ROOT / 'output/teacher_review_20260906.zip')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--verify', type=Path)
    args = parser.parse_args()
    if args.verify:
        result = verify(args.verify.resolve())
    else:
        root = args.root.resolve()
        entries, missing = collect(root)
        verify_thermal_ledger(root, entries)
        result = dict(mode='plan', source_files=len(entries),
                      source_bytes=sum(e['size'] for e in entries), missing=missing,
                      selected_directories=list(DIRECTORIES))
        if not args.plan:
            if missing:
                raise FileNotFoundError('缺少必需材料：' + ', '.join(missing))
            result = package(root, args.output.resolve(), entries)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
