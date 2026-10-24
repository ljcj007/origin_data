#!/usr/bin/env python3
"""生成或核对结果目录的 SHA-256 账本。"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest().upper()


def collect(root: Path, output: Path) -> list[dict[str, int | str]]:
    files = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.resolve() == output.resolve():
            continue
        files.append(
            {
                "path": path.relative_to(root).as_posix(),
                "size": path.stat().st_size,
                "sha256": digest(path),
            }
        )
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"结果目录不存在：{root}")
    output = (args.output or root / "checksums.sha256.json").resolve()
    entries = collect(root, output)

    if args.verify:
        if not output.is_file():
            raise SystemExit(f"校验账本不存在：{output}")
        expected = json.loads(output.read_text(encoding="utf-8"))
        if expected.get("files") != entries:
            raise SystemExit(f"校验失败：{root}")
        print(f"校验通过：{root}，{len(entries)} 个文件")
        return

    payload = {
        "schema_version": 1,
        "root": root.name,
        "file_count": len(entries),
        "files": entries,
    }
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"账本已写入：{output}，{len(entries)} 个文件")


if __name__ == "__main__":
    main()
