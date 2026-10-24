#!/usr/bin/env python3
"""汇总 SSE 有限尺寸扫描并检查两个独立副本。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SIZES = (4, 6, 8, 10, 12, 14, 16)
EXPECTED_REPLICAS = (1, 2)
EXPECTED_CELL_COUNT = len(EXPECTED_SIZES) * len(EXPECTED_REPLICAS)
EXPECTED_BETA_OVER_L2 = 2.0
EXPECTED_BINS = 64
EXPECTED_EQUIL_SWEEPS = 10000
EXPECTED_BIN_SWEEPS = 1000
FIT_KEYS = (
    "energy_per_site",
    "staggered_m2",
    "ferrimagnetic_m2",
    "sublattice_plus_m2",
    "sublattice_minus_m2",
    "sublattice_cross",
    "spin_stiffness_cell_average",
)
REQUIRED_BIN_COLUMNS = (
    "energy_per_site",
    "staggered_m2",
    "staggered_m4",
    "ferrimagnetic_m2",
    "sublattice_plus_m2",
    "sublattice_minus_m2",
    "sublattice_cross",
    "spin_stiffness_cell_x",
    "spin_stiffness_cell_y",
    "spin_stiffness_cell_average",
    "spin_stiffness_times_L",
    "winding_squared",
    "mean_expansion_order",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def mean_error(values: list[float]) -> tuple[float, float]:
    return statistics.fmean(values), statistics.stdev(values) / math.sqrt(len(values))


def lag1(values: list[float]) -> float:
    mean = statistics.fmean(values)
    denominator = sum((value - mean) ** 2 for value in values)
    if denominator == 0.0:
        return 0.0
    numerator = sum(
        (values[index] - mean) * (values[index - 1] - mean)
        for index in range(1, len(values))
    )
    return numerator / denominator


def binder_jackknife(m2: list[float], m4: list[float]) -> tuple[float, float]:
    count = len(m2)

    def ratio(mean_m2: float, mean_m4: float) -> float:
        return 1.0 - 3.0 * mean_m4 / mean_m2**2

    value = ratio(statistics.fmean(m2), statistics.fmean(m4))
    total_m2, total_m4 = sum(m2), sum(m4)
    samples = [
        ratio((total_m2 - m2[index]) / (count - 1), (total_m4 - m4[index]) / (count - 1))
        for index in range(count)
    ]
    center = statistics.fmean(samples)
    error = math.sqrt((count - 1) / count * sum((sample - center) ** 2 for sample in samples))
    return value, error


def load_bins(path: Path, expected_count: int) -> dict[str, list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    if len(rows) != expected_count:
        raise RuntimeError(f"bin 数量应为 {expected_count}，实际为 {len(rows)}：{path}")
    if len(fieldnames) != len(set(fieldnames)):
        raise RuntimeError(f"bin 表存在重复列名：{path}")
    missing = [name for name in ("bin", *REQUIRED_BIN_COLUMNS) if name not in fieldnames]
    if missing:
        raise RuntimeError(f"bin 表缺少列 {missing}：{path}")
    try:
        indices = [int(row["bin"]) for row in rows]
        columns = {
            name: [float(row[name]) for row in rows]
            for name in fieldnames
            if name != "bin"
        }
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"bin 表含非数值字段：{path}") from exc
    if indices != list(range(1, expected_count + 1)):
        raise RuntimeError(f"bin 编号不连续：{path}")
    if any(not math.isfinite(value) for values in columns.values() for value in values):
        raise RuntimeError(f"bin 表含非有限值：{path}")
    return columns


def same_value(actual: Any, expected: Any) -> bool:
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        try:
            return math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=1.0e-12)
        except (TypeError, ValueError):
            return False
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            same_value(item, target) for item, target in zip(actual, expected)
        )
    return actual == expected


def blocked_error(replica_values: list[list[float]], block_size: int) -> float:
    blocks: list[float] = []
    for values in replica_values:
        if len(values) % block_size:
            raise RuntimeError(f"bin 数量不能被块长 {block_size} 整除")
        blocks.extend(
            statistics.fmean(values[start : start + block_size])
            for start in range(0, len(values), block_size)
        )
    return statistics.stdev(blocks) / math.sqrt(len(blocks))


def conservative_error(values: list[float]) -> float:
    """计算单副本保守分块误差。"""
    errors = [mean_error(values)[1]]
    errors.extend(blocked_error([values], block) for block in (2, 4, 8))
    return max(errors)


def parse_input(path: Path) -> dict[str, Any]:
    """读取正式 SSE 输入。"""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) != 5:
        raise RuntimeError(f"SSE 输入行数异常：{path}")
    try:
        lx, ly = (int(value) for value in lines[1].split())
        beta, type_count = lines[2].split()
        couplings = [float(value) for value in lines[3].split()]
        nbins, equil, bin_sweeps, seed = (int(value) for value in lines[4].split())
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"SSE 输入格式异常：{path}") from exc
    return {
        "graph": lines[0].replace("\\", "/"),
        "Lx": lx,
        "Ly": ly,
        "beta": float(beta),
        "type_count": int(type_count),
        "couplings": couplings,
        "nbins": nbins,
        "equil_sweeps": equil,
        "bin_sweeps": bin_sweeps,
        "seed": seed,
    }


def validate_protocol(run_dir: Path, spec: dict[str, Any]) -> None:
    """锁定正式 14 单元协议和本地来源。"""
    cells = spec.get("cells", [])
    ids = [str(cell.get("cell_id")) for cell in cells]
    if len(cells) != EXPECTED_CELL_COUNT or len(set(ids)) != EXPECTED_CELL_COUNT:
        raise RuntimeError("正式协议必须包含 14 个唯一单元")
    observed = {
        (int(cell["params"]["L"]), int(cell["params"]["replica"]))
        for cell in cells
    }
    expected = {(size, replica) for size in EXPECTED_SIZES for replica in EXPECTED_REPLICAS}
    if observed != expected:
        raise RuntimeError("正式尺寸或副本轴不完整")
    if any(
        not same_value(cell["params"].get("beta_over_L2"), EXPECTED_BETA_OVER_L2)
        for cell in cells
    ):
        raise RuntimeError("正式 beta J/L^2 规则不一致")

    settings = spec.get("settings", {})
    expected_settings = {
        "model": "spin_half_antiferromagnetic_heisenberg",
        "graph": "graphs/super_honeycomb_poscar.dat",
        "boundary": "periodic_xy",
        "couplings": [1.0],
        "basis_sites": 5,
        "basis_plus": 2,
        "basis_minus": 3,
        "nbins": EXPECTED_BINS,
        "equil_sweeps": EXPECTED_EQUIL_SWEEPS,
        "bin_sweeps": EXPECTED_BIN_SWEEPS,
    }
    for key, expected_value in expected_settings.items():
        if not same_value(settings.get(key), expected_value):
            raise RuntimeError(f"正式运行设置不一致：{key}")

    provenance = spec.get("provenance", {})
    config_pairs = (
        (run_dir / "config" / "settings.json", settings, "settings"),
        (run_dir / "config" / "provenance.json", provenance, "provenance"),
    )
    for path, expected_value, label in config_pairs:
        if not path.is_file() or json.loads(path.read_text(encoding="utf-8")) != expected_value:
            raise RuntimeError(f"config/{label}.json 与 run spec 不一致")

    local_sources = {
        "poscar_sha256": ROOT / "POSCAR (2).vasp",
        "graph_sha256": ROOT / str(settings["graph"]),
        "source_sha256": ROOT / "sse_graph_stiffness_v2.f90",
    }
    for key, path in local_sources.items():
        if sha256_file(path) != str(provenance.get(key, "")).upper():
            raise RuntimeError(f"本地来源哈希不一致：{key}")


def validate_input(
    input_data: dict[str, Any], manifest: dict[str, Any], spec: dict[str, Any], cell_id: str
) -> None:
    """核对输入文件与清单参数。"""
    params = manifest["params"]
    settings = spec["settings"]
    size = int(params["L"])
    expected = {
        "Lx": size,
        "Ly": size,
        "beta": float(params["beta"]),
        "type_count": len(settings["couplings"]),
        "couplings": settings["couplings"],
        "nbins": int(settings["nbins"]),
        "equil_sweeps": int(settings["equil_sweeps"]),
        "bin_sweeps": int(settings["bin_sweeps"]),
        "seed": int(manifest["execution"]["seed"]),
    }
    for key, expected_value in expected.items():
        if not same_value(input_data.get(key), expected_value):
            raise RuntimeError(f"输入字段 {key} 与清单不一致：{cell_id}")
    graph_suffix = str(settings["graph"]).replace("\\", "/")
    if not str(input_data["graph"]).endswith(graph_suffix):
        raise RuntimeError(f"输入格点图路径不一致：{cell_id}")


def validate_bin_identities(
    bins: dict[str, list[float]], spec: dict[str, Any], size: int, beta: float, cell_id: str
) -> None:
    """核对逐 bin 的能量、磁矩和刚度恒等式。"""
    n_plus = int(spec["settings"]["basis_plus"]) * size * size
    n_minus = int(spec["settings"]["basis_minus"]) * size * size
    nsite = n_plus + n_minus
    coupling_sum = 6.0 * size * size
    tolerance = 5.0e-12
    for index in range(len(bins["energy_per_site"])):
        energy = -bins["mean_expansion_order"][index] / (beta * nsite) + coupling_sum / (4.0 * nsite)
        plus_weight = (n_plus / nsite) ** 2 * bins["sublattice_plus_m2"][index]
        minus_weight = (n_minus / nsite) ** 2 * bins["sublattice_minus_m2"][index]
        cross_weight = 2.0 * n_plus * n_minus / nsite**2 * bins["sublattice_cross"][index]
        expected = {
            "energy_per_site": energy,
            "ferrimagnetic_m2": plus_weight + minus_weight + cross_weight,
            "staggered_m2": plus_weight + minus_weight - cross_weight,
            "spin_stiffness_cell_average": 0.5
            * (bins["spin_stiffness_cell_x"][index] + bins["spin_stiffness_cell_y"][index]),
            "spin_stiffness_times_L": bins["spin_stiffness_cell_average"][index] * size,
            "winding_relation": 0.75 * bins["winding_squared"][index] / beta,
        }
        actual = {
            name: bins[name][index]
            for name in (
                "energy_per_site",
                "ferrimagnetic_m2",
                "staggered_m2",
                "spin_stiffness_cell_average",
                "spin_stiffness_times_L",
            )
        }
        actual["winding_relation"] = bins["spin_stiffness_cell_average"][index]
        for name, expected_value in expected.items():
            if not math.isclose(actual[name], expected_value, rel_tol=0.0, abs_tol=tolerance):
                raise RuntimeError(f"逐-bin 恒等式失败 {name}：{cell_id}/bin-{index + 1}")


def validate_cell(run_dir: Path, cell: dict[str, Any], spec: dict[str, Any]) -> dict[str, Any]:
    cell_dir = run_dir / "cells" / cell["cell_id"]
    manifest_path = cell_dir / "manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"缺少单元：{cell['cell_id']}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise RuntimeError(f"不支持的清单版本：{cell['cell_id']}")
    if not manifest.get("status", {}).get("success"):
        raise RuntimeError(f"失败单元：{cell['cell_id']}")
    if manifest.get("cell_id") != cell["cell_id"]:
        raise RuntimeError(f"单元编号不一致：{cell['cell_id']}")
    manifest_settings = manifest.get("settings", {})
    for key, expected in spec["settings"].items():
        if not same_value(manifest_settings.get(key), expected):
            raise RuntimeError(f"运行设置 {key} 不一致：{cell['cell_id']}")
    for key, expected in cell.get("params", {}).items():
        actual = manifest.get("params", {}).get(key)
        if not same_value(actual, expected):
            raise RuntimeError(f"参数 {key} 不一致：{cell['cell_id']}")
    params = manifest.get("params", {})
    expected_beta = float(params["beta_over_L2"]) * int(params["L"]) ** 2
    if not same_value(params.get("beta"), expected_beta):
        raise RuntimeError(f"beta != beta_over_L2*L^2：{cell['cell_id']}")
    provenance = spec["provenance"]
    manifest_provenance = manifest.get("provenance", {})
    for key, expected in provenance.items():
        if not same_value(manifest_provenance.get(key), expected):
            raise RuntimeError(f"清单 {key} 不一致：{cell['cell_id']}")
    execution = manifest.get("execution", {})
    if execution.get("returncode") != 0:
        raise RuntimeError(f"运行返回码非零：{cell['cell_id']}")
    if not isinstance(execution.get("seed"), int):
        raise RuntimeError(f"缺少整数随机种子：{cell['cell_id']}")
    for key in ("graph_sha256", "executable_sha256"):
        if execution.get(key, "").upper() != provenance[key].upper():
            raise RuntimeError(f"{key} 不一致：{cell['cell_id']}")
    artifact_hashes = manifest.get("artifact_sha256", {})
    for name, filename in (
        ("input", "sse_graph.in"),
        ("bins", "sse_graph_bins.csv"),
        ("summary", "sse_graph_summary.txt"),
    ):
        path = cell_dir / filename
        if sha256_file(path) != artifact_hashes.get(name, "").upper():
            raise RuntimeError(f"{name} 工件哈希不一致：{cell['cell_id']}")
    validate_input(parse_input(cell_dir / "sse_graph.in"), manifest, spec, cell["cell_id"])
    size = int(params["L"])
    basis_sites = int(spec["settings"]["basis_sites"])
    imbalance = abs(int(spec["settings"]["basis_minus"]) - int(spec["settings"]["basis_plus"]))
    expected_spin = 0.5 * imbalance * size * size
    expected_target = expected_spin * (expected_spin + 1.0) / (basis_sites * size * size) ** 2
    benchmarks = manifest.get("benchmarks", {})
    if not same_value(benchmarks.get("lieb_mattis_ground_spin"), expected_spin):
        raise RuntimeError(f"Lieb-Mattis 总自旋不一致：{cell['cell_id']}")
    if not same_value(benchmarks.get("zero_temperature_ferrimagnetic_m2"), expected_target):
        raise RuntimeError(f"Lieb-Mattis 磁矩目标不一致：{cell['cell_id']}")
    return manifest


def replica_z(replica_bins: list[dict[str, list[float]]], key: str) -> float:
    means = [statistics.fmean(bins[key]) for bins in replica_bins]
    errors = [conservative_error(bins[key]) for bins in replica_bins]
    sigma = math.hypot(*errors)
    return abs(means[0] - means[1]) / sigma if sigma else math.inf


def summarize(run_dir: Path) -> tuple[list[dict[str, float | int]], dict[str, Any]]:
    spec = json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))
    validate_protocol(run_dir, spec)
    batch_path = run_dir / "direct_batch_status.json"
    if not batch_path.is_file():
        raise RuntimeError("缺少批处理完成状态")
    batch = json.loads(batch_path.read_text(encoding="utf-8"))
    if batch.get("run_id") != spec.get("run_id"):
        raise RuntimeError("批处理 run_id 不一致")
    if batch.get("workers") != len(spec["cells"]):
        raise RuntimeError("批处理并发数与正式协议不一致")
    batch_cells = batch.get("cells", [])
    expected_ids = {cell["cell_id"] for cell in spec["cells"]}
    if len(batch_cells) != EXPECTED_CELL_COUNT or {
        cell.get("cell_id") for cell in batch_cells
    } != expected_ids:
        raise RuntimeError("批处理单元集合不完整")
    if any(cell.get("returncode") != 0 for cell in batch_cells):
        raise RuntimeError("批处理中存在失败单元")
    groups: dict[int, list[tuple[dict[str, Any], dict[str, list[float]]]]] = {}
    for cell in spec["cells"]:
        manifest = validate_cell(run_dir, cell, spec)
        size = int(manifest["params"]["L"])
        bins = load_bins(
            run_dir / "cells" / cell["cell_id"] / "sse_graph_bins.csv",
            int(spec["settings"]["nbins"]),
        )
        validate_bin_identities(
            bins,
            spec,
            int(manifest["params"]["L"]),
            float(manifest["params"]["beta"]),
            cell["cell_id"],
        )
        for key in (*FIT_KEYS, "staggered_m4"):
            value, error = mean_error(bins[key])
            recorded = manifest["observables"][key]
            if not same_value(recorded["value"], value) or not same_value(recorded["error"], error):
                raise RuntimeError(f"观测量 {key} 与原始 bin 不一致：{cell['cell_id']}")
        binder, binder_error = binder_jackknife(bins["staggered_m2"], bins["staggered_m4"])
        recorded_binder = manifest["observables"]["binder_z_staggered"]
        if not same_value(recorded_binder["value"], binder) or not same_value(
            recorded_binder["error"], binder_error
        ):
            raise RuntimeError(f"Binder 与原始 bin 不一致：{cell['cell_id']}")
        for diagnostic, key in (
            ("energy_bin_lag1", "energy_per_site"),
            ("staggered_m2_bin_lag1", "staggered_m2"),
            ("ferrimagnetic_m2_bin_lag1", "ferrimagnetic_m2"),
        ):
            if not same_value(manifest["diagnostics"][diagnostic], lag1(bins[key])):
                raise RuntimeError(f"诊断量 {diagnostic} 与原始 bin 不一致：{cell['cell_id']}")
        groups.setdefault(size, []).append((manifest, bins))

    rows: list[dict[str, float | int]] = []
    replica_checks: list[dict[str, float | int]] = []
    replica_lag1: list[dict[str, float | str]] = []
    for size in sorted(groups):
        group = groups[size]
        if len(group) != 2:
            raise RuntimeError(f"L={size} 的副本数不是 2")
        manifests = [item[0] for item in group]
        if len({int(manifest["execution"]["seed"]) for manifest in manifests}) != 2:
            raise RuntimeError(f"L={size} 的两个副本使用了相同种子")
        if not same_value(manifests[0]["params"]["beta"], manifests[1]["params"]["beta"]):
            raise RuntimeError(f"L={size} 的两个副本 beta 不一致")
        pooled = {
            name: [value for _, bins in group for value in bins[name]]
            for name in group[0][1]
        }
        params = manifests[0]["params"]
        row: dict[str, float | int] = {
            "L": size,
            "N": int(spec["settings"]["basis_sites"]) * size * size,
            "beta": float(params["beta"]),
            "replicas": len(group),
            "pooled_bins": len(pooled["energy_per_site"]),
        }
        for key in FIT_KEYS:
            value, naive_error = mean_error(pooled[key])
            block_errors = {
                block: blocked_error([item[1][key] for item in group], block)
                for block in (2, 4, 8)
            }
            error = max(naive_error, *block_errors.values())
            row[key] = value
            row[f"{key}_error"] = error
            row[f"{key}_error_naive"] = naive_error
            for block, block_error in block_errors.items():
                row[f"{key}_error_block{block}"] = block_error
            row[f"{key}_blocking_ratio_max"] = error / naive_error if naive_error else math.inf
        binder, binder_error = binder_jackknife(pooled["staggered_m2"], pooled["staggered_m4"])
        row["binder_z_staggered"] = binder
        row["binder_z_staggered_error"] = binder_error
        spin = 0.5 * abs(
            int(spec["settings"]["basis_minus"]) - int(spec["settings"]["basis_plus"])
        ) * size * size
        target = spin * (spin + 1.0) / (int(spec["settings"]["basis_sites"]) * size * size) ** 2
        row["lieb_mattis_target"] = target
        row["ferrimagnetic_ratio"] = float(row["ferrimagnetic_m2"]) / target
        row["ferrimagnetic_ratio_error"] = float(row["ferrimagnetic_m2_error"]) / target
        for column, key in (
            ("energy_bin_lag1_max_replica", "energy_per_site"),
            ("staggered_m2_bin_lag1_max_replica", "staggered_m2"),
            ("ferrimagnetic_m2_bin_lag1_max_replica", "ferrimagnetic_m2"),
        ):
            replica_values = [lag1(item[1][key]) for item in group]
            row[column] = max(replica_values, key=abs)
        rows.append(row)
        replica_checks.append(
            {
                "L": size,
                **{f"{key}_z": replica_z([item[1] for item in group], key) for key in FIT_KEYS},
            }
        )
        for manifest in manifests:
            replica_lag1.append(
                {
                    "cell_id": manifest["cell_id"],
                    "L": size,
                    **{
                        key: float(manifest["diagnostics"][key])
                        for key in (
                            "energy_bin_lag1",
                            "staggered_m2_bin_lag1",
                            "ferrimagnetic_m2_bin_lag1",
                        )
                    },
                }
            )

    z_values = [
        float(value)
        for check in replica_checks
        for key, value in check.items()
        if key.endswith("_z")
    ]
    all_seeds = [
        int(manifest["execution"]["seed"])
        for group in groups.values()
        for manifest, _ in group
    ]
    if len(set(all_seeds)) != len(all_seeds):
        raise RuntimeError("正式扫描存在重复随机种子")
    summary = {
        "run_id": spec["run_id"],
        "batch_workers": int(batch["workers"]),
        "batch_elapsed_seconds": float(batch["elapsed_seconds"]),
        "cell_count": len(spec["cells"]),
        "size_count": len(rows),
        "sizes": [int(row["L"]) for row in rows],
        "replicas_per_size": 2,
        "bins_per_size": [int(row["pooled_bins"]) for row in rows],
        "max_abs_recomputed_replica_lag1": max(
            abs(float(row[key]))
            for row in rows
            for key in (
                "energy_bin_lag1_max_replica",
                "staggered_m2_bin_lag1_max_replica",
                "ferrimagnetic_m2_bin_lag1_max_replica",
            )
        ),
        "max_abs_replica_lag1": max(
            abs(float(value))
            for item in replica_lag1
            for key, value in item.items()
            if key.endswith("_lag1")
        ),
        "max_replica_difference_sigma": max(z_values),
        "replica_difference_error_policy": "each replica uses max(naive, block2, block4, block8) SEM",
        "max_blocking_error_ratio": max(
            float(value)
            for row in rows
            for key, value in row.items()
            if key.endswith("_blocking_ratio_max")
        ),
        "fit_error_policy": "max(naive, block2, block4, block8) replica-aware standard error",
        "replica_checks": replica_checks,
        "replica_lag1": replica_lag1,
    }
    return rows, summary


def write_csv(rows: list[dict[str, float | int]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    rows, summary = summarize(run_dir)
    write_csv(rows, run_dir / "fss_observables.csv")
    (run_dir / "fss_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
