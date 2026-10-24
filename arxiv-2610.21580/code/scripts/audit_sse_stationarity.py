#!/usr/bin/env python3
"""独立审计 SSE bin 的平稳性与副本一致性。"""

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
EXPECTED_CELL_COUNT = 14
EXPECTED_BIN_COUNT = 64
BLOCK_SIZES = (1, 2, 4, 8)
SPLIT_BLOCK_SIZE = 4
ALPHA = 0.05
OBSERVABLES = (
    "energy_per_site",
    "staggered_m2",
    "ferrimagnetic_m2",
    "sublattice_plus_m2",
    "sublattice_minus_m2",
    "sublattice_cross",
    "spin_stiffness_cell_average",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def block_sem(values: list[float], block_size: int) -> float:
    """计算连续分块均值的标准误。"""
    if block_size < 1 or len(values) % block_size:
        raise ValueError("分块长度必须整除样本数")
    blocks = [
        statistics.fmean(values[start : start + block_size])
        for start in range(0, len(values), block_size)
    ]
    if len(blocks) < 2:
        raise ValueError("分块后样本数不足")
    return statistics.stdev(blocks) / math.sqrt(len(blocks))


def sem_profile(values: list[float]) -> dict[str, float]:
    """汇总多种分块标准误。"""
    errors = {f"block{size}": block_sem(values, size) for size in BLOCK_SIZES}
    conservative = max(errors.values())
    naive = errors["block1"]
    return {
        **errors,
        "conservative": conservative,
        "blocking_ratio": conservative / naive if naive else math.inf,
    }


def lag1(values: list[float]) -> float:
    """计算相邻 bin 的一阶相关系数。"""
    mean = statistics.fmean(values)
    denominator = sum((value - mean) ** 2 for value in values)
    if denominator == 0.0:
        return 0.0
    numerator = sum(
        (values[index] - mean) * (values[index - 1] - mean)
        for index in range(1, len(values))
    )
    return numerator / denominator


def signed_z(delta: float, first_error: float, second_error: float) -> float:
    denominator = math.hypot(first_error, second_error)
    if denominator:
        return delta / denominator
    if delta == 0.0:
        return 0.0
    return math.copysign(math.inf, delta)


def two_sided_normal_p(z_value: float) -> float:
    if math.isinf(z_value):
        return 0.0
    return math.erfc(abs(z_value) / math.sqrt(2.0))


def split_half_statistics(values: list[float]) -> dict[str, float]:
    """比较前后两个等长区段。"""
    if len(values) < 16 or len(values) % 2:
        raise ValueError("前后半段分析需要偶数且不少于 16 个 bin")
    midpoint = len(values) // 2
    first = values[:midpoint]
    second = values[midpoint:]
    first_profile = sem_profile(first)
    second_profile = sem_profile(second)
    first_mean = statistics.fmean(first)
    second_mean = statistics.fmean(second)
    delta = second_mean - first_mean
    z_block4 = signed_z(
        delta,
        first_profile[f"block{SPLIT_BLOCK_SIZE}"],
        second_profile[f"block{SPLIT_BLOCK_SIZE}"],
    )
    z_conservative = signed_z(
        delta,
        first_profile["conservative"],
        second_profile["conservative"],
    )
    return {
        "first_mean": first_mean,
        "second_mean": second_mean,
        "delta": delta,
        "first_sem_block4": first_profile["block4"],
        "second_sem_block4": second_profile["block4"],
        "z_block4": z_block4,
        "p_block4": two_sided_normal_p(z_block4),
        "first_sem_conservative": first_profile["conservative"],
        "second_sem_conservative": second_profile["conservative"],
        "z_conservative": z_conservative,
        "p_conservative": two_sided_normal_p(z_conservative),
    }


def read_bins(path: Path, expected_count: int = EXPECTED_BIN_COUNT) -> dict[str, list[float]]:
    """读取并严格检查顺序 bin。"""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if len(rows) != expected_count:
        raise ValueError(f"bin 数量应为 {expected_count}，实际为 {len(rows)}：{path}")
    if [int(row["bin"]) for row in rows] != list(range(1, expected_count + 1)):
        raise ValueError(f"bin 编号不连续：{path}")
    missing = [name for name in OBSERVABLES if name not in fieldnames]
    if missing:
        raise ValueError(f"bin 表缺少列 {missing}：{path}")
    columns = {name: [float(row[name]) for row in rows] for name in OBSERVABLES}
    if any(not math.isfinite(value) for values in columns.values() for value in values):
        raise ValueError(f"bin 表含非有限值：{path}")
    return columns


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def discover_completed_cells(run_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """按 run spec 发现已完整下载的成功单元。"""
    spec = load_json(run_dir / "run_spec.json")
    cells = spec.get("cells", [])
    if len(cells) != EXPECTED_CELL_COUNT:
        raise ValueError("run spec 不是预期的 14 单元 SSE 扫描")
    expected_pairs = {
        (size, replica) for size in (4, 6, 8, 10, 12, 14, 16) for replica in (1, 2)
    }
    actual_pairs = {
        (int(cell["params"]["L"]), int(cell["params"]["replica"])) for cell in cells
    }
    if actual_pairs != expected_pairs:
        raise ValueError("run spec 的尺寸或副本轴异常")

    completed: list[dict[str, Any]] = []
    for planned in cells:
        cell_id = str(planned["cell_id"])
        cell_dir = run_dir / "cells" / cell_id
        manifest_path = cell_dir / "manifest.json"
        bins_path = cell_dir / "sse_graph_bins.csv"
        if not manifest_path.is_file() and not bins_path.is_file():
            continue
        if not manifest_path.is_file() or not bins_path.is_file():
            raise ValueError(f"单元文件不完整：{cell_id}")
        manifest = load_json(manifest_path)
        if manifest.get("cell_id") != cell_id:
            raise ValueError(f"单元编号不一致：{cell_id}")
        if not manifest.get("status", {}).get("success"):
            raise ValueError(f"单元未成功：{cell_id}")
        params = manifest.get("params", {})
        if int(params.get("L")) != int(planned["params"]["L"]):
            raise ValueError(f"单元尺寸不一致：{cell_id}")
        if int(params.get("replica")) != int(planned["params"]["replica"]):
            raise ValueError(f"单元副本不一致：{cell_id}")
        expected_hash = str(manifest.get("artifact_sha256", {}).get("bins", "")).upper()
        if not expected_hash or sha256_file(bins_path) != expected_hash:
            raise ValueError(f"bin 文件哈希不一致：{cell_id}")
        completed.append(
            {
                "cell_id": cell_id,
                "L": int(params["L"]),
                "replica": int(params["replica"]),
                "seed": int(manifest["execution"]["seed"]),
                "bins_path": bins_path,
                "bins_sha256": expected_hash,
            }
        )
    completed.sort(key=lambda item: (item["L"], item["replica"]))
    return spec, completed


def audit(run_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """执行独立诊断，不产生排除决定。"""
    spec, cells = discover_completed_cells(run_dir)
    bins_by_cell: dict[str, dict[str, list[float]]] = {}
    cell_rows: list[dict[str, Any]] = []
    for cell in cells:
        bins = read_bins(cell["bins_path"])
        bins_by_cell[cell["cell_id"]] = bins
        for observable in OBSERVABLES:
            values = bins[observable]
            profile = sem_profile(values)
            split = split_half_statistics(values)
            cell_rows.append(
                {
                    "cell_id": cell["cell_id"],
                    "L": cell["L"],
                    "replica": cell["replica"],
                    "observable": observable,
                    "n_bins": len(values),
                    "mean": statistics.fmean(values),
                    "sem_naive": profile["block1"],
                    "sem_block2": profile["block2"],
                    "sem_block4": profile["block4"],
                    "sem_block8": profile["block8"],
                    "sem_conservative": profile["conservative"],
                    "blocking_ratio": profile["blocking_ratio"],
                    "lag1": lag1(values),
                    **{f"split_{key}": value for key, value in split.items()},
                }
            )

    groups: dict[int, list[dict[str, Any]]] = {}
    for cell in cells:
        groups.setdefault(int(cell["L"]), []).append(cell)
    replica_rows: list[dict[str, Any]] = []
    for size, group in sorted(groups.items()):
        if len(group) != 2:
            continue
        group.sort(key=lambda item: item["replica"])
        if [item["replica"] for item in group] != [1, 2]:
            raise ValueError(f"L={size} 的副本编号异常")
        if group[0]["seed"] == group[1]["seed"]:
            raise ValueError(f"L={size} 的两个副本种子相同")
        for observable in OBSERVABLES:
            first_values = bins_by_cell[group[0]["cell_id"]][observable]
            second_values = bins_by_cell[group[1]["cell_id"]][observable]
            first_profile = sem_profile(first_values)
            second_profile = sem_profile(second_values)
            first_mean = statistics.fmean(first_values)
            second_mean = statistics.fmean(second_values)
            delta = second_mean - first_mean
            z_block4 = signed_z(
                delta, first_profile["block4"], second_profile["block4"]
            )
            z_conservative = signed_z(
                delta, first_profile["conservative"], second_profile["conservative"]
            )
            replica_rows.append(
                {
                    "L": size,
                    "observable": observable,
                    "replica1_cell": group[0]["cell_id"],
                    "replica2_cell": group[1]["cell_id"],
                    "replica1_mean": first_mean,
                    "replica2_mean": second_mean,
                    "delta_replica2_minus_replica1": delta,
                    "replica1_sem_block4": first_profile["block4"],
                    "replica2_sem_block4": second_profile["block4"],
                    "z_block4": z_block4,
                    "p_block4": two_sided_normal_p(z_block4),
                    "replica1_sem_conservative": first_profile["conservative"],
                    "replica2_sem_conservative": second_profile["conservative"],
                    "z_conservative": z_conservative,
                    "p_conservative": two_sided_normal_p(z_conservative),
                }
            )

    split_test_count = len(cell_rows)
    replica_test_count = len(replica_rows)
    for row in cell_rows:
        row["nominal_split_attention"] = abs(float(row["split_z_block4"])) >= 3.0
        row["familywise_split_attention"] = (
            float(row["split_p_block4"]) <= ALPHA / split_test_count
            if split_test_count
            else False
        )
        row["lag1_attention"] = abs(float(row["lag1"])) >= 0.5
        row["blocking_attention"] = float(row["blocking_ratio"]) >= 2.0
    for row in replica_rows:
        row["nominal_replica_attention"] = abs(float(row["z_conservative"])) >= 3.0
        row["familywise_replica_attention"] = (
            float(row["p_conservative"]) <= ALPHA / replica_test_count
            if replica_test_count
            else False
        )

    expected_ids = [str(cell["cell_id"]) for cell in spec["cells"]]
    available_ids = [str(cell["cell_id"]) for cell in cells]
    missing_ids = [cell_id for cell_id in expected_ids if cell_id not in set(available_ids)]
    complete_sizes = sorted(size for size, group in groups.items() if len(group) == 2)
    max_split_block4 = max(cell_rows, key=lambda row: abs(float(row["split_z_block4"])))
    max_split_conservative = max(
        cell_rows, key=lambda row: abs(float(row["split_z_conservative"]))
    )
    max_replica = max(replica_rows, key=lambda row: abs(float(row["z_conservative"])))
    max_lag1 = max(cell_rows, key=lambda row: abs(float(row["lag1"])))
    max_blocking = max(cell_rows, key=lambda row: float(row["blocking_ratio"]))
    summary = {
        "schema_version": 1,
        "audit_id": "sse_fss_stationarity_v1",
        "source_run_id": spec.get("run_id"),
        "status": "complete_review" if len(cells) == EXPECTED_CELL_COUNT else "provisional_review",
        "provisional": len(cells) != EXPECTED_CELL_COUNT,
        "scope": "independent read-only diagnostic; not a preregistered acceptance gate",
        "decision_policy": "report only; never exclude or rerun data automatically",
        "completeness": {
            "expected_cell_count": EXPECTED_CELL_COUNT,
            "available_cell_count": len(cells),
            "available_cell_ids": available_ids,
            "missing_cell_ids": missing_ids,
            "complete_sizes": complete_sizes,
        },
        "methods": {
            "observables": list(OBSERVABLES),
            "bins_per_cell": EXPECTED_BIN_COUNT,
            "split": "bins 1-32 versus 33-64",
            "split_primary_error": "4-bin contiguous blocking SEM within each half",
            "conservative_error": "max SEM from contiguous block sizes 1, 2, 4, 8",
            "lag1_definition": "adjacent-bin covariance divided by full-bin sum of squares",
            "nominal_attention": "absolute z >= 3; descriptive only",
            "familywise_attention": "two-sided normal p <= 0.05/number_of_tests; descriptive only",
            "lag1_attention": "absolute lag1 >= 0.5; descriptive only",
            "blocking_attention": "conservative SEM / naive SEM >= 2; descriptive only",
        },
        "test_counts": {
            "split_half": split_test_count,
            "replica": replica_test_count,
        },
        "attention_counts": {
            "nominal_split": sum(bool(row["nominal_split_attention"]) for row in cell_rows),
            "familywise_split": sum(
                bool(row["familywise_split_attention"]) for row in cell_rows
            ),
            "nominal_replica": sum(
                bool(row["nominal_replica_attention"]) for row in replica_rows
            ),
            "familywise_replica": sum(
                bool(row["familywise_replica_attention"]) for row in replica_rows
            ),
            "lag1": sum(bool(row["lag1_attention"]) for row in cell_rows),
            "blocking": sum(bool(row["blocking_attention"]) for row in cell_rows),
        },
        "extrema": {
            "max_abs_split_z_block4": select_fields(
                max_split_block4, "cell_id", "L", "replica", "observable", "split_z_block4"
            ),
            "max_abs_split_z_conservative": select_fields(
                max_split_conservative,
                "cell_id",
                "L",
                "replica",
                "observable",
                "split_z_conservative",
            ),
            "max_abs_replica_z_conservative": select_fields(
                max_replica, "L", "observable", "z_conservative"
            ),
            "max_abs_lag1": select_fields(
                max_lag1, "cell_id", "L", "replica", "observable", "lag1"
            ),
            "max_blocking_ratio": select_fields(
                max_blocking,
                "cell_id",
                "L",
                "replica",
                "observable",
                "blocking_ratio",
            ),
        },
        "source_artifacts": [
            {
                "cell_id": cell["cell_id"],
                "L": cell["L"],
                "replica": cell["replica"],
                "bins_sha256": cell["bins_sha256"],
            }
            for cell in cells
        ],
        "audit_source_sha256": sha256_file(Path(__file__)),
    }
    return cell_rows, replica_rows, summary


def select_fields(row: dict[str, Any], *names: str) -> dict[str, Any]:
    return {name: row[name] for name in names}


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError("没有可写入的诊断行")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def format_extreme(item: dict[str, Any], value_name: str) -> str:
    labels = []
    if "cell_id" in item:
        labels.append(str(item["cell_id"]))
    if "L" in item:
        labels.append(f"L={item['L']}")
    if "replica" in item:
        labels.append(f"副本 {item['replica']}")
    if "observable" in item:
        labels.append(str(item["observable"]))
    return f"{', '.join(labels)}：{float(item[value_name]):.6g}"


def render_report(summary: dict[str, Any]) -> str:
    status = "临时报告（provisional）" if summary["provisional"] else "完整报告"
    completeness = summary["completeness"]
    extrema = summary["extrema"]
    attention = summary["attention_counts"]
    missing = ", ".join(completeness["missing_cell_ids"]) or "无"
    return f"""# SSE 有限尺寸数据平稳性审计 v1

状态：**{status}**。已读取 {completeness['available_cell_count']}/{completeness['expected_cell_count']} 个单元；缺失：{missing}。

## 边界

这是独立、只读的探索性诊断，不属于预注册正式分析门。脚本只报告异常线索，不自动排除数据、不自动要求重跑，也不改动 `results/sse_fss_v1`。若完整数据出现需关注项，应结合时间序列、多重比较和物理量相关性人工判断。

## 方法

- 对 7 个正式拟合观测量逐副本比较 bins 1–32 与 33–64。
- 固定 4-bin 连续分块标准误给出主要前后半段 z 值；另以 1、2、4、8-bin 标准误最大值给出保守 z 值。
- 对同一 L 的两个独立副本，用各自保守标准误计算一致性 z 值。
- 同时报告 lag-1 相关和“保守标准误/朴素标准误”阻塞比。
- |z|≥3、|lag-1|≥0.5、阻塞比≥2 只作提示；另列 0.05 族错误率的 Bonferroni 提示。

## 当前极值

- 固定 4-bin 前后半段：{format_extreme(extrema['max_abs_split_z_block4'], 'split_z_block4')}
- 保守前后半段：{format_extreme(extrema['max_abs_split_z_conservative'], 'split_z_conservative')}
- 保守副本差：{format_extreme(extrema['max_abs_replica_z_conservative'], 'z_conservative')}
- 最大 |lag-1|：{format_extreme(extrema['max_abs_lag1'], 'lag1')}
- 最大阻塞比：{format_extreme(extrema['max_blocking_ratio'], 'blocking_ratio')}

## 提示计数

- 前后半段名义提示：{attention['nominal_split']}；Bonferroni 族内提示：{attention['familywise_split']}
- 副本差名义提示：{attention['nominal_replica']}；Bonferroni 族内提示：{attention['familywise_replica']}
- lag-1 提示：{attention['lag1']}；阻塞比提示：{attention['blocking']}

逐项数值见 `stationarity_cells.csv` 与 `replica_consistency.csv`。14/14 单元到齐后使用相同命令重跑即可，报告状态会自动切换为 `complete_review`。
"""


def write_checksum_ledger(output_dir: Path) -> None:
    records = []
    ledger_path = output_dir / "checksums.sha256.json"
    for path in sorted(item for item in output_dir.iterdir() if item.is_file()):
        if path == ledger_path:
            continue
        records.append(
            {
                "path": path.name,
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    payload = {
        "schema_version": 1,
        "root": output_dir.name,
        "file_count": len(records),
        "files": records,
    }
    ledger_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_outputs(
    output_dir: Path,
    cell_rows: list[dict[str, Any]],
    replica_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(cell_rows, output_dir / "stationarity_cells.csv")
    write_csv(replica_rows, output_dir / "replica_consistency.csv")
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(render_report(summary), encoding="utf-8")
    write_checksum_ledger(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "results" / "sse_fss_v1")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results" / "sse_fss_stationarity_v1",
    )
    args = parser.parse_args()
    cell_rows, replica_rows, summary = audit(args.run_dir.resolve())
    write_outputs(args.output_dir.resolve(), cell_rows, replica_rows, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
