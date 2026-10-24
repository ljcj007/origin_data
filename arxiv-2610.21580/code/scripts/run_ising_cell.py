#!/usr/bin/env python3
"""执行一个 Ising 参数扫描单元并生成清单。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dump_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def deterministic_seed(run_id: str, cell_id: str) -> int:
    digest = hashlib.sha256(f"{run_id}:{cell_id}".encode()).digest()
    return int.from_bytes(digest[:4], "big") % 2_000_000_000 + 1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def mean_error(values: list[float]) -> tuple[float, float]:
    mean = statistics.fmean(values)
    error = statistics.stdev(values) / math.sqrt(len(values))
    return mean, error


def lag1(values: list[float]) -> float:
    mean = statistics.fmean(values)
    variance = sum((value - mean) ** 2 for value in values)
    if variance == 0.0:
        return 0.0
    return sum((values[i] - mean) * (values[i - 1] - mean) for i in range(1, len(values))) / variance


def binder_jackknife(q2: list[float], q4: list[float]) -> tuple[float, float]:
    count = len(q2)
    mean_q2 = statistics.fmean(q2)
    mean_q4 = statistics.fmean(q4)
    value = 1.0 - mean_q4 / (3.0 * mean_q2**2)
    samples = []
    sum_q2, sum_q4 = sum(q2), sum(q4)
    for index in range(count):
        leave_q2 = (sum_q2 - q2[index]) / (count - 1)
        leave_q4 = (sum_q4 - q4[index]) / (count - 1)
        samples.append(1.0 - leave_q4 / (3.0 * leave_q2**2))
    center = statistics.fmean(samples)
    error = math.sqrt((count - 1) / count * sum((sample - center) ** 2 for sample in samples))
    return value, error


def parse_bins(path: Path) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < 8:
        raise ValueError("独立 bin 数量不足")
    columns = {name: [float(row[name]) for row in rows] for name in rows[0] if name != "bin"}
    observables = {}
    for name, values in columns.items():
        value, error = mean_error(values)
        observables[name] = {"value": value, "error": error}
    binder, binder_error = binder_jackknife(columns["staggered_m2"], columns["staggered_m4"])
    observables["binder_staggered"] = {
        "value": binder,
        "error": binder_error,
        "estimator": "delete-one-bin jackknife",
    }
    diagnostics = {
        "energy_bin_lag1": lag1(columns["energy_per_site"]),
        "staggered_m2_bin_lag1": lag1(columns["staggered_m2"]),
    }
    return observables, diagnostics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-spec", type=Path, required=True)
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--executable", type=Path, default=ROOT / "ising_graph.exe")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    spec = load_json(args.run_spec.resolve())
    cells = {cell["cell_id"]: cell for cell in spec["cells"]}
    if args.cell_id not in cells:
        raise SystemExit("run spec 中没有该 cell")
    cell = cells[args.cell_id]
    params = cell.get("params", {})
    settings = {**spec.get("settings", {}), **cell.get("settings", {})}
    run_dir = ROOT / spec["run_dir"].replace("\\", "/")
    cell_dir = run_dir / "cells" / args.cell_id
    cell_dir.mkdir(parents=True, exist_ok=True)

    size = int(params["L"])
    temperature = float(params["temperature"])
    budget = int(settings["budget"])
    nbins = int(settings["nbins"])
    if budget % nbins != 0:
        raise SystemExit("budget 必须能被 nbins 整除")
    seed = deterministic_seed(spec["run_id"], args.cell_id)
    graph = (ROOT / settings["graph"]).resolve()
    executable = args.executable.resolve()
    actual_graph_hash = sha256_file(graph)
    actual_executable_hash = sha256_file(executable)
    provenance = spec.get("provenance", {})
    if actual_graph_hash != provenance.get("graph_sha256", "").upper():
        raise SystemExit("格点图哈希与运行清单不一致")
    if actual_executable_hash != provenance.get("executable_sha256", "").upper():
        raise SystemExit("可执行文件哈希与运行清单不一致")
    input_path = cell_dir / "ising_graph.in"
    input_path.write_text(
        f"{graph}\n{size} {size}\n{temperature:.17g} {settings['field']:.17g} {settings['coupling']:.17g}\n"
        f"{nbins} {settings['equil_sweeps']} {budget // nbins} {seed}\n"
        f"{settings.get('algorithm', 'metropolis')}\n",
        encoding="utf-8",
    )

    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable), str(input_path.resolve())],
        cwd=cell_dir,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    (cell_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (cell_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    success = completed.returncode == 0 and (cell_dir / "ising_bins.csv").is_file()
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "cell_id": args.cell_id,
        "params": params,
        "settings": settings,
        "provenance": provenance,
        "status": {"success": success},
        "execution": {
            "returncode": completed.returncode,
            "seed": seed,
            "elapsed_seconds": elapsed,
            "graph_sha256": actual_graph_hash,
            "executable_sha256": actual_executable_hash,
        },
        "artifacts": {
            "bins": str((cell_dir / "ising_bins.csv").resolve()),
            "summary": str((cell_dir / "ising_summary.txt").resolve()),
        },
    }
    if success:
        observables, diagnostics = parse_bins(cell_dir / "ising_bins.csv")
        manifest["observables"] = observables
        manifest["diagnostics"] = diagnostics
        manifest["artifact_sha256"] = {
            "input": sha256_file(input_path),
            "bins": sha256_file(cell_dir / "ising_bins.csv"),
            "summary": sha256_file(cell_dir / "ising_summary.txt"),
        }
    dump_json(cell_dir / "manifest.json", manifest)
    if not success:
        raise SystemExit(f"cell 失败，返回码 {completed.returncode}")
    print(cell_dir / "manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
