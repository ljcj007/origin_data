#!/usr/bin/env python3
"""合并动力学 SSE 副本并审计虚时相关函数。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np


def read_modes(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_tgrid(path: Path) -> np.ndarray:
    values = path.read_text(encoding="utf-8").split()
    count = int(values[0]) + 1
    tau = np.asarray([float(value) for value in values[1:]], dtype=float)
    if len(tau) != count:
        raise ValueError("虚时网格长度错误")
    return tau


def read_cor(path: Path, nmode: int, ntau: int) -> tuple[np.ndarray, list[tuple[int, int]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    block = ntau + 1
    if len(lines) % (nmode * block):
        raise ValueError("动力学相关函数块数错误")
    nbin = len(lines) // (nmode * block)
    values = np.empty((nbin, nmode, ntau), dtype=float)
    headers: list[tuple[int, int]] = []
    cursor = 0
    for bin_index in range(nbin):
        for mode_index in range(nmode):
            header = tuple(int(value) for value in lines[cursor].split())
            cursor += 1
            if len(header) != 2:
                raise ValueError("动力学相关函数头错误")
            if bin_index == 0:
                headers.append(header)
            elif header != headers[mode_index]:
                raise ValueError("副本内部模式顺序不一致")
            values[bin_index, mode_index] = [
                float(value) for value in lines[cursor : cursor + ntau]
            ]
            cursor += ntau
    return values, headers


def standard_error(values: np.ndarray) -> np.ndarray:
    if values.shape[0] < 2:
        return np.full(values.shape[1:], math.nan)
    return np.std(values, axis=0, ddof=1) / math.sqrt(values.shape[0])


def write_cor(path: Path, arrays: list[np.ndarray], headers: list[tuple[int, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for values in arrays:
            for one_bin in values:
                for mode_index, (mode_id, q_label) in enumerate(headers):
                    handle.write(f"{mode_id:8d}{q_label:8d}\n")
                    for value in one_bin[mode_index]:
                        handle.write(f"{value:26.16E}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("replicas", nargs="+", type=Path)
    args = parser.parse_args()

    replicas = [path.resolve() for path in args.replicas]
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"输出目录已存在：{output}")
    if len(replicas) < 2:
        raise ValueError("至少需要两个独立副本")

    mode_bytes = (replicas[0] / "dynamics_modes.csv").read_bytes()
    tgrid_bytes = (replicas[0] / "dynamics_tgrid.dat").read_bytes()
    modes = read_modes(replicas[0] / "dynamics_modes.csv")
    tau = read_tgrid(replicas[0] / "dynamics_tgrid.dat")
    arrays: list[np.ndarray] = []
    headers: list[tuple[int, int]] | None = None
    for replica in replicas:
        if (replica / "dynamics_modes.csv").read_bytes() != mode_bytes:
            raise ValueError("副本模式清单不一致")
        if (replica / "dynamics_tgrid.dat").read_bytes() != tgrid_bytes:
            raise ValueError("副本虚时网格不一致")
        values, current_headers = read_cor(
            replica / "dynamics_cor.dat", len(modes), len(tau)
        )
        if headers is None:
            headers = current_headers
        elif current_headers != headers:
            raise ValueError("副本模式标签不一致")
        arrays.append(values)

    normalized = [values / values[:, :, :1] for values in arrays]
    pair_records = []
    for left_index in range(len(normalized)):
        for right_index in range(left_index + 1, len(normalized)):
            left = normalized[left_index]
            right = normalized[right_index]
            left_mean = np.mean(left, axis=0)
            right_mean = np.mean(right, axis=0)
            left_error = standard_error(left)
            right_error = standard_error(right)
            denominator = np.sqrt(left_error**2 + right_error**2)
            difference = np.abs(left_mean - right_mean)
            zscore = np.divide(
                difference,
                denominator,
                out=np.zeros_like(difference),
                where=denominator > 0,
            )
            for mode_index, row in enumerate(modes):
                valid = (
                    (left_mean[mode_index] > 0)
                    & (right_mean[mode_index] > 0)
                    & (left_error[mode_index] / left_mean[mode_index] <= 0.1)
                    & (right_error[mode_index] / right_mean[mode_index] <= 0.1)
                )
                invalid = np.flatnonzero(~valid[1:])
                stop = int(invalid[0] + 1) if invalid.size else len(tau)
                usable = zscore[mode_index, 1:stop]
                if usable.size:
                    maximum_index = int(np.argmax(usable)) + 1
                    maximum_z = float(zscore[mode_index, maximum_index])
                    points_over = int(np.count_nonzero(usable > 3.0))
                else:
                    maximum_index = 0
                    maximum_z = 0.0
                    points_over = 0
                pair_records.append(
                    {
                        "replica_left": left_index + 1,
                        "replica_right": right_index + 1,
                        "mode_id": int(row["mode_id"]),
                        "q_index": int(row["q_index"]),
                        "channel": row["channel"],
                        "max_abs_z": maximum_z,
                        "tau_at_max": float(tau[maximum_index]),
                        "usable_tau_max": float(tau[max(stop - 1, 0)]),
                        "points_over_3sigma": points_over,
                    }
                )

    beta_values = []
    for replica in replicas:
        summary = (replica / "sse_graph_summary.txt").read_text(encoding="utf-8")
        match = re.search(r"beta\s*=\s*([0-9Ee+\-.]+)", summary)
        if not match:
            raise ValueError("无法读取 beta")
        beta_values.append(float(match.group(1)))
    if not np.allclose(beta_values, beta_values[0], rtol=0.0, atol=1.0e-12):
        raise ValueError("副本 beta 不一致")

    output.mkdir(parents=True)
    shutil.copyfile(replicas[0] / "dynamics_modes.csv", output / "dynamics_modes.csv")
    shutil.copyfile(replicas[0] / "dynamics_tgrid.dat", output / "dynamics_tgrid.dat")
    write_cor(output / "dynamics_cor.dat", arrays, headers or [])
    (output / "sse_graph_summary.txt").write_text(
        f"pooled dynamics replicas = {len(replicas)}\n"
        f"beta = {beta_values[0]:.17g}\n"
        f"bins = {sum(values.shape[0] for values in arrays)}\n",
        encoding="utf-8",
    )
    audit = {
        "format": "sac-replica-audit-v1",
        "replicas": [str(path) for path in replicas],
        "bins_per_replica": [int(values.shape[0]) for values in arrays],
        "pooled_bins": int(sum(values.shape[0] for values in arrays)),
        "beta": beta_values[0],
        "max_abs_z": max(record["max_abs_z"] for record in pair_records),
        "modes_with_any_point_over_3sigma": sum(
            record["points_over_3sigma"] > 0 for record in pair_records
        ),
        "records": pair_records,
    }
    (output / "replica_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: audit[key] for key in ("pooled_bins", "max_abs_z", "modes_with_any_point_over_3sigma")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
