#!/usr/bin/env python3
"""在直连服务器上并行执行尚未完成的扫描单元。"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def cell_succeeded(run_dir: Path, cell_id: str) -> bool:
    path = run_dir / "cells" / cell_id / "manifest.json"
    if not path.is_file():
        return False
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return bool(manifest.get("status", {}).get("success"))


def run_cell(
    run_spec: Path, executable: Path, cell_runner: Path, cell_id: str
) -> tuple[str, int, str, float]:
    command = [
        sys.executable,
        str(cell_runner),
        "--run-spec",
        str(run_spec),
        "--cell-id",
        cell_id,
        "--executable",
        str(executable),
    ]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - started
    message = completed.stdout.strip() or completed.stderr.strip()
    return cell_id, completed.returncode, message, elapsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-spec", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--cell-runner", type=Path, default=ROOT / "scripts" / "run_ising_cell.py")
    parser.add_argument("--workers", type=int, default=10)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    run_spec = args.run_spec.resolve()
    executable = args.executable.resolve()
    cell_runner = args.cell_runner.resolve()
    spec = json.loads(run_spec.read_text(encoding="utf-8"))
    run_dir = ROOT / spec["run_dir"].replace("\\", "/")
    completed = [cell["cell_id"] for cell in spec["cells"] if cell_succeeded(run_dir, cell["cell_id"])]
    pending = [cell["cell_id"] for cell in spec["cells"] if cell["cell_id"] not in completed]
    started = time.time()
    results = [
        (cell_id, 0, str(run_dir / "cells" / cell_id / "manifest.json"), 0.0)
        for cell_id in completed
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [
            executor.submit(run_cell, run_spec, executable, cell_runner, cell_id)
            for cell_id in pending
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result[0]} returncode={result[1]}", flush=True)
    status = {
        "run_id": spec["run_id"],
        "workers": args.workers,
        "elapsed_seconds": time.time() - started,
        "cells": [
            {
                "cell_id": cell_id,
                "returncode": returncode,
                "message": message,
                "elapsed_seconds": elapsed,
            }
            for cell_id, returncode, message, elapsed in sorted(results)
        ],
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "direct_batch_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return 0 if all(returncode == 0 for _, returncode, _, _ in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
