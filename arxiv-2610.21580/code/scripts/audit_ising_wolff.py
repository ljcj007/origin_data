#!/usr/bin/env python3
"""严格复核正式 Wolff Ising 扫描。"""

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
BLOCK_SIZES = (1, 2, 4, 8)
EXPECTED_SIZES = (8, 12, 16, 24, 32, 48, 64)
EXPECTED_TEMPERATURES = (0.996, 1.000, 1.004, 1.008, 1.012)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def close(a: float, b: float, tolerance: float = 1.0e-12) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def lag1(values: list[float]) -> float:
    mean = statistics.fmean(values)
    variance = sum((value - mean) ** 2 for value in values)
    if variance == 0.0:
        return 0.0
    return sum(
        (values[index] - mean) * (values[index - 1] - mean)
        for index in range(1, len(values))
    ) / variance


def mean_error(values: list[float], block_size: int = 1) -> tuple[float, float]:
    if block_size < 1 or len(values) % block_size:
        raise ValueError("分块大小必须整除样本数")
    blocks = [
        statistics.fmean(values[start : start + block_size])
        for start in range(0, len(values), block_size)
    ]
    if len(blocks) < 2:
        raise ValueError("分块后样本数不足")
    return statistics.fmean(values), statistics.stdev(blocks) / math.sqrt(len(blocks))


def binder_jackknife(
    q2: list[float], q4: list[float], block_size: int = 1
) -> tuple[float, float]:
    if len(q2) != len(q4) or block_size < 1 or len(q2) % block_size:
        raise ValueError("Binder 分块输入不合法")
    block_q2 = [
        statistics.fmean(q2[start : start + block_size])
        for start in range(0, len(q2), block_size)
    ]
    block_q4 = [
        statistics.fmean(q4[start : start + block_size])
        for start in range(0, len(q4), block_size)
    ]
    count = len(block_q2)
    if count < 2:
        raise ValueError("Binder 分块后样本数不足")
    mean_q2 = statistics.fmean(q2)
    mean_q4 = statistics.fmean(q4)
    value = 1.0 - mean_q4 / (3.0 * mean_q2**2)
    sum_q2, sum_q4 = sum(block_q2), sum(block_q4)
    samples = []
    for index in range(count):
        leave_q2 = (sum_q2 - block_q2[index]) / (count - 1)
        leave_q4 = (sum_q4 - block_q4[index]) / (count - 1)
        samples.append(1.0 - leave_q4 / (3.0 * leave_q2**2))
    center = statistics.fmean(samples)
    error = math.sqrt(
        (count - 1) / count * sum((sample - center) ** 2 for sample in samples)
    )
    return value, error


def split_half_z(values: list[float], block_size: int = 4) -> float:
    midpoint = len(values) // 2
    first, second = values[:midpoint], values[midpoint:]
    first_mean, first_error = mean_error(first, block_size)
    second_mean, second_error = mean_error(second, block_size)
    denominator = math.hypot(first_error, second_error)
    difference = second_mean - first_mean
    if denominator == 0.0:
        return 0.0 if difference == 0.0 else math.copysign(math.inf, difference)
    return difference / denominator


def binder_split_half_z(q2: list[float], q4: list[float], block_size: int = 4) -> float:
    midpoint = len(q2) // 2
    first, first_error = binder_jackknife(q2[:midpoint], q4[:midpoint], block_size)
    second, second_error = binder_jackknife(q2[midpoint:], q4[midpoint:], block_size)
    denominator = math.hypot(first_error, second_error)
    difference = second - first
    if denominator == 0.0:
        return 0.0 if difference == 0.0 else math.copysign(math.inf, difference)
    return difference / denominator


def read_bins(path: Path) -> dict[str, list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 64 or [int(row["bin"]) for row in rows] != list(range(1, 65)):
        raise ValueError(f"{path} 的 64 个 bin 不完整")
    columns = {
        name: [float(row[name]) for row in rows]
        for name in rows[0]
        if name != "bin"
    }
    if any(not math.isfinite(value) for values in columns.values() for value in values):
        raise ValueError(f"{path} 含非有限值")
    return columns


def crossing_roots(
    rows: list[dict[str, float | int | str]], error_key: str
) -> list[dict[str, Any]]:
    sizes = sorted({int(row["L"]) for row in rows})
    crossings: list[dict[str, Any]] = []
    for small, large in zip(sizes, sizes[1:]):
        first = {float(row["temperature"]): row for row in rows if int(row["L"]) == small}
        second = {float(row["temperature"]): row for row in rows if int(row["L"]) == large}
        differences = []
        for temperature in sorted(set(first) & set(second)):
            delta = float(second[temperature]["binder"]) - float(first[temperature]["binder"])
            sigma = math.hypot(
                float(second[temperature][error_key]), float(first[temperature][error_key])
            )
            differences.append((temperature, delta, sigma))
        roots = []
        for left, right in zip(differences, differences[1:]):
            t1, d1, s1 = left
            t2, d2, s2 = right
            if d1 * d2 > 0.0 or d2 == d1:
                continue
            denominator = d2 - d1
            root = (t1 * d2 - t2 * d1) / denominator
            derivative1 = d2 * (t1 - t2) / denominator**2
            derivative2 = d1 * (t2 - t1) / denominator**2
            roots.append(
                {
                    "temperature": root,
                    "statistical_error": math.hypot(derivative1 * s1, derivative2 * s2),
                    "bracket": [t1, t2],
                }
            )
        crossings.append({"L_small": small, "L_large": large, "roots": roots})
    return crossings


def verify_checksum_ledger(run_dir: Path) -> int:
    ledger = load_json(run_dir / "checksums.sha256.json")
    records = ledger["files"]
    if int(ledger["file_count"]) != len(records):
        raise ValueError("校验账本文件数不一致")
    seen: set[str] = set()
    for record in records:
        relative = str(record["path"])
        if relative in seen:
            raise ValueError(f"校验账本重复路径：{relative}")
        seen.add(relative)
        path = run_dir / relative
        if not path.is_file():
            raise ValueError(f"校验账本缺失文件：{relative}")
        if path.stat().st_size != int(record["size"]):
            raise ValueError(f"文件大小不一致：{relative}")
        if sha256_file(path) != str(record["sha256"]).upper():
            raise ValueError(f"文件哈希不一致：{relative}")
    actual = {
        path.relative_to(run_dir).as_posix()
        for path in run_dir.rglob("*")
        if path.is_file() and path.name != "checksums.sha256.json"
    }
    if seen != actual:
        raise ValueError("校验账本未精确覆盖结果目录")
    return len(records)


def parse_input(path: Path) -> dict[str, Any]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) != 5:
        raise ValueError(f"输入文件格式异常：{path}")
    lx, ly = (int(value) for value in lines[1].split())
    temperature, field, coupling = (float(value) for value in lines[2].split())
    nbins, equil, bin_sweeps, seed = (int(value) for value in lines[3].split())
    return {
        "graph": lines[0],
        "Lx": lx,
        "Ly": ly,
        "temperature": temperature,
        "field": field,
        "coupling": coupling,
        "nbins": nbins,
        "equil_sweeps": equil,
        "bin_sweeps": bin_sweeps,
        "seed": seed,
        "algorithm": lines[4],
    }


def audit(run_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    spec = load_json(run_dir / "run_spec.json")
    settings = spec["settings"]
    provenance = spec["provenance"]
    expected_pairs = {
        (temperature, size)
        for temperature in EXPECTED_TEMPERATURES
        for size in EXPECTED_SIZES
    }
    observed_pairs = {
        (float(cell["params"]["temperature"]), int(cell["params"]["L"]))
        for cell in spec["cells"]
    }
    if observed_pairs != expected_pairs or len(spec["cells"]) != 35:
        raise ValueError("正式温度或尺寸网格不完整")
    if settings != load_json(run_dir / "config" / "settings.json"):
        raise ValueError("settings 配置与 run spec 不一致")
    if provenance != load_json(run_dir / "config" / "provenance.json"):
        raise ValueError("provenance 配置与 run spec 不一致")

    local_hashes = {
        "poscar_sha256": sha256_file(ROOT / "POSCAR (2).vasp"),
        "graph_sha256": sha256_file(ROOT / str(settings["graph"])),
        "source_sha256": sha256_file(ROOT / "ising_graph.f90"),
    }
    for name, value in local_hashes.items():
        if str(provenance[name]).upper() != value:
            raise ValueError(f"本地 {name} 与正式运行不一致")

    cell_dirs = {path.name for path in (run_dir / "cells").iterdir() if path.is_dir()}
    expected_ids = {str(cell["cell_id"]) for cell in spec["cells"]}
    if cell_dirs != expected_ids:
        raise ValueError("cell 目录集合与 run spec 不一致")

    rows: list[dict[str, Any]] = []
    seeds: set[int] = set()
    execution_hash_count = 0
    for cell in spec["cells"]:
        cell_id = str(cell["cell_id"])
        cell_dir = run_dir / "cells" / cell_id
        manifest = load_json(cell_dir / "manifest.json")
        if manifest["cell_id"] != cell_id or manifest["params"] != cell["params"]:
            raise ValueError(f"cell 参数不一致：{cell_id}")
        if manifest["settings"] != settings or manifest["provenance"] != provenance:
            raise ValueError(f"cell 配置或来源不一致：{cell_id}")
        if not manifest["status"]["success"] or int(manifest["execution"]["returncode"]) != 0:
            raise ValueError(f"cell 未成功：{cell_id}")
        if "executable_sha256" in manifest["execution"]:
            execution_hash_count += 1

        actual_input = parse_input(cell_dir / "ising_graph.in")
        size = int(cell["params"]["L"])
        temperature = float(cell["params"]["temperature"])
        expected_input = {
            "Lx": size,
            "Ly": size,
            "temperature": temperature,
            "field": float(settings["field"]),
            "coupling": float(settings["coupling"]),
            "nbins": int(settings["nbins"]),
            "equil_sweeps": int(settings["equil_sweeps"]),
            "bin_sweeps": int(settings["budget"]) // int(settings["nbins"]),
            "algorithm": str(settings["algorithm"]),
        }
        for name, expected in expected_input.items():
            actual = actual_input[name]
            if isinstance(expected, float):
                if not close(float(actual), expected):
                    raise ValueError(f"输入字段不一致：{cell_id}/{name}")
            elif actual != expected:
                raise ValueError(f"输入字段不一致：{cell_id}/{name}")
        seed = int(actual_input["seed"])
        if seed != int(manifest["execution"]["seed"]) or seed in seeds:
            raise ValueError(f"随机种子重复或不一致：{cell_id}")
        seeds.add(seed)

        bins = read_bins(cell_dir / "ising_bins.csv")
        if any(not close(value, 1.0) for value in bins["clusters_per_sweep"]):
            raise ValueError(f"并非每步固定一个 Wolff 集群：{cell_id}")
        q2, q4 = bins["staggered_m2"], bins["staggered_m4"]
        binder, original_error = binder_jackknife(q2, q4)
        recorded = manifest["observables"]["binder_staggered"]
        if not close(binder, float(recorded["value"])) or not close(
            original_error, float(recorded["error"])
        ):
            raise ValueError(f"Binder 无法由原始 bin 复现：{cell_id}")
        energy, energy_error = mean_error(bins["energy_per_site"])
        recorded_energy = manifest["observables"]["energy_per_site"]
        if not close(energy, float(recorded_energy["value"])) or not close(
            energy_error, float(recorded_energy["error"])
        ):
            raise ValueError(f"能量无法由原始 bin 复现：{cell_id}")

        binder_errors = {
            block_size: binder_jackknife(q2, q4, block_size)[1]
            for block_size in BLOCK_SIZES
        }
        energy_errors = {
            block_size: mean_error(bins["energy_per_site"], block_size)[1]
            for block_size in BLOCK_SIZES
        }
        conservative_binder = max(binder_errors.values())
        conservative_energy = max(energy_errors.values())
        measured_lag = lag1(q2)
        recorded_lag = float(manifest["diagnostics"]["staggered_m2_bin_lag1"])
        if not close(measured_lag, recorded_lag):
            raise ValueError(f"lag-1 无法由原始 bin 复现：{cell_id}")
        rows.append(
            {
                "cell_id": cell_id,
                "temperature": temperature,
                "L": size,
                "binder": binder,
                "binder_error_original": original_error,
                "binder_error_block2": binder_errors[2],
                "binder_error_block4": binder_errors[4],
                "binder_error_block8": binder_errors[8],
                "binder_error_conservative": conservative_binder,
                "binder_block_ratio": conservative_binder / original_error,
                "energy": energy,
                "energy_error_original": energy_error,
                "energy_error_conservative": conservative_energy,
                "energy_block_ratio": conservative_energy / energy_error,
                "lag1_staggered_m2": measured_lag,
                "binder_first_second_half_z": binder_split_half_z(q2, q4),
                "energy_first_second_half_z": split_half_z(bins["energy_per_site"]),
            }
        )

    scan_rows = {}
    with (run_dir / "scan_observables.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            scan_rows[str(row["cell_id"])] = row
    if set(scan_rows) != expected_ids:
        raise ValueError("scan_observables.csv 单元集合不一致")
    for row in rows:
        recorded = scan_rows[row["cell_id"]]
        if not close(float(recorded["binder"]), float(row["binder"])) or not close(
            float(recorded["binder_error"]), float(row["binder_error_original"])
        ):
            raise ValueError(f"扫描汇总与原始 bin 不一致：{row['cell_id']}")

    conservative_crossings = crossing_roots(rows, "binder_error_conservative")
    exact_temperature = 2.0 / math.acosh(2.0 + math.sqrt(3.0))
    largest = next(
        item for item in conservative_crossings
        if item["L_small"] == 48 and item["L_large"] == 64
    )["roots"]
    if len(largest) != 1:
        raise ValueError("最大尺寸 crossing 根数量异常")
    largest_root = largest[0]
    largest_root["difference_from_exact"] = largest_root["temperature"] - exact_temperature
    largest_root["statistical_z"] = (
        largest_root["difference_from_exact"] / largest_root["statistical_error"]
    )

    max_binder_row = max(rows, key=lambda row: float(row["binder_block_ratio"]))
    max_energy_row = max(rows, key=lambda row: float(row["energy_block_ratio"]))
    max_binder_half_row = max(
        rows, key=lambda row: abs(float(row["binder_first_second_half_z"]))
    )
    max_energy_half_row = max(
        rows, key=lambda row: abs(float(row["energy_first_second_half_z"]))
    )
    report = {
        "schema_version": 1,
        "status": "pass_with_provenance_boundary",
        "run_id": spec["run_id"],
        "cell_count": len(rows),
        "bin_count_per_cell": 64,
        "measurements_per_bin": int(settings["budget"]) // int(settings["nbins"]),
        "total_measurements": len(rows) * int(settings["budget"]),
        "unique_seed_count": len(seeds),
        "checksum_ledger_file_count": verify_checksum_ledger(run_dir),
        "local_hashes_verified": local_hashes,
        "manifest_execution_hash_count": execution_hash_count,
        "executable_provenance_boundary": (
            "35 个 manifest 继承同一个已声明的 Linux 可执行文件哈希，"
            "但 execution 字段没有记录运行时重新计算的哈希；本地没有该 Linux 二进制，"
            "因此当前只能确认声明一致，不能独立重算二进制哈希。"
        ),
        "exact_temperature": exact_temperature,
        "conservative_crossings": conservative_crossings,
        "largest_pair_crossing": largest_root,
        "max_abs_staggered_m2_lag1": max(abs(float(row["lag1_staggered_m2"])) for row in rows),
        "max_binder_block_ratio": {
            "value": max_binder_row["binder_block_ratio"],
            "cell_id": max_binder_row["cell_id"],
            "L": max_binder_row["L"],
            "temperature": max_binder_row["temperature"],
        },
        "max_energy_block_ratio": {
            "value": max_energy_row["energy_block_ratio"],
            "cell_id": max_energy_row["cell_id"],
            "L": max_energy_row["L"],
            "temperature": max_energy_row["temperature"],
        },
        "stationarity_diagnostic": {
            "definition": "前后各 32 个 bin，分别按 4-bin 分块后的均值差 z 值；仅作描述性检查",
            "max_abs_binder_half_z": {
                "value": abs(float(max_binder_half_row["binder_first_second_half_z"])),
                "cell_id": max_binder_half_row["cell_id"],
            },
            "max_abs_energy_half_z": {
                "value": abs(float(max_energy_half_row["energy_first_second_half_z"])),
                "cell_id": max_energy_half_row["cell_id"],
            },
            "independent_hot_cold_replica_available": False,
        },
        "uncertainty_boundary": {
            "statistical_error": "相邻 1/2/4/8 个原始 bin 分块所得标准误差的最大值",
            "temperature_grid_step": 0.004,
            "crossing_interpolation": "相邻温度点线性插值",
            "temperature_grid_systematic_included": False,
            "finite_size_systematic_included": False,
            "precision_extrapolation_claim_allowed": False,
            "allowed_interpretation": "对精确临界温度和格点实现的数值验证",
        },
    }
    return report, sorted(rows, key=lambda row: (float(row["temperature"]), int(row["L"])))


def write_outputs(output_dir: Path, report: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "audit_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    with (output_dir / "audit_observables.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir", type=Path, default=ROOT / "results" / "ising_tc_wolff_v1"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results" / "ising_tc_wolff_audit_v1"
    )
    args = parser.parse_args()
    report, rows = audit(args.run_dir.resolve())
    write_outputs(args.output_dir.resolve(), report, rows)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
