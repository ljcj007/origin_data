#!/usr/bin/env python3
"""验证正式 SSE 下载暂存区，避免失败工件进入正式目录。"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ARTIFACTS = {
    "input": "sse_graph.in",
    "bins": "sse_graph_bins.csv",
    "summary": "sse_graph_summary.txt",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate(run_spec: Path, stage: Path, target_ids: list[str]) -> dict[str, Any]:
    spec = load_json(run_spec)
    batch = load_json(stage / "direct_batch_status.json")
    expected = {cell["cell_id"]: cell for cell in spec["cells"]}
    expected_ids = sorted(expected)
    batch_ids = sorted(str(cell["cell_id"]) for cell in batch["cells"])
    if batch.get("run_id") != spec.get("run_id"):
        raise ValueError("批次 run_id 与登记不符")
    if batch_ids != expected_ids or len(batch_ids) != 14:
        raise ValueError("批次没有恰好覆盖登记的 14 个单元")
    if any(int(cell["returncode"]) != 0 for cell in batch["cells"]):
        raise ValueError("批次状态包含失败单元")
    if sorted(target_ids) != ["cell-0013", "cell-0014"]:
        raise ValueError("暂存验收仅允许两个登记的 L=16 单元")
    if not (stage / "direct_batch.log").is_file():
        raise FileNotFoundError("暂存区缺少 direct_batch.log")

    verified = []
    for cell_id in target_ids:
        registered = expected[cell_id]
        cell_dir = stage / cell_id
        manifest = load_json(cell_dir / "manifest.json")
        if manifest.get("cell_id") != cell_id:
            raise ValueError(f"单元编号不符：{cell_id}")
        if manifest.get("status", {}).get("success") is not True:
            raise ValueError(f"单元未成功：{cell_id}")
        if int(manifest.get("execution", {}).get("returncode", -1)) != 0:
            raise ValueError(f"单元返回码失败：{cell_id}")
        for key in ("L", "replica", "beta_over_L2"):
            if manifest["params"].get(key) != registered["params"].get(key):
                raise ValueError(f"单元参数不符：{cell_id}/{key}")
        expected_beta = float(registered["params"]["beta_over_L2"]) * int(
            registered["params"]["L"]
        ) ** 2
        if float(manifest["params"].get("beta", -1.0)) != expected_beta:
            raise ValueError(f"单元 beta 不符：{cell_id}")
        if manifest.get("settings") != spec.get("settings"):
            raise ValueError(f"单元设置不符：{cell_id}")
        if manifest.get("provenance") != spec.get("provenance"):
            raise ValueError(f"单元来源不符：{cell_id}")
        execution = manifest["execution"]
        provenance = spec["provenance"]
        if execution.get("graph_sha256") != provenance.get("graph_sha256"):
            raise ValueError(f"运行图哈希不符：{cell_id}")
        if execution.get("executable_sha256") != provenance.get("executable_sha256"):
            raise ValueError(f"运行程序哈希不符：{cell_id}")
        for key, name in ARTIFACTS.items():
            path = cell_dir / name
            if not path.is_file():
                raise FileNotFoundError(f"暂存工件缺失：{cell_id}/{name}")
            expected_hash = str(manifest["artifact_sha256"].get(key, "")).upper()
            if sha256_file(path) != expected_hash:
                raise ValueError(f"暂存工件哈希不符：{cell_id}/{name}")
        verified.append(cell_id)

    return {
        "status": "passed",
        "run_id": spec["run_id"],
        "batch_cells": len(batch_ids),
        "verified_staged_cells": verified,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-spec", type=Path, required=True)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--cell-id", action="append", dest="cell_ids", required=True)
    args = parser.parse_args()
    result = validate(args.run_spec.resolve(), args.stage.resolve(), args.cell_ids)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
