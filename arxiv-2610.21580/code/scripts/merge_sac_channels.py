#!/usr/bin/env python3
"""合并同一 QMC 数据的分通道 SAC 结果。"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()

    summaries = []
    spectra = []
    seen_records: set[tuple[int, str]] = set()
    seen_spectra: set[tuple[int, str, float]] = set()
    records = []
    for directory in args.inputs:
        directory = directory.resolve()
        summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
        summaries.append(summary)
        for record in summary["records"]:
            key = (int(record["mode_id"]), str(record["channel"]))
            if key in seen_records:
                raise ValueError(f"重复模式：{key}")
            seen_records.add(key)
            records.append(record)
        with (directory / "spectra.csv").open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                key = (int(row["mode_id"]), row["channel"], float(row["omega"]))
                if key in seen_spectra:
                    raise ValueError(f"重复谱点：{key}")
                seen_spectra.add(key)
                spectra.append(row)

    args.output.mkdir(parents=True, exist_ok=False)
    payload = {
        "format": "sac-batch-merged-v1",
        "sources": [str(path.resolve()) for path in args.inputs],
        "beta": summaries[0]["beta"],
        "settings": summaries[0]["settings"],
        "counts": {
            "selected": len(records),
            "complete": sum(record["status"] == "complete" for record in records),
            "skipped_or_failed": sum(record["status"] != "complete" for record in records),
        },
        "records": sorted(records, key=lambda record: int(record["mode_id"])),
    }
    (args.output / "summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    spectra.sort(key=lambda row: (int(row["mode_id"]), float(row["omega"])))
    with (args.output / "spectra.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(spectra[0]))
        writer.writeheader()
        writer.writerows(spectra)
    print(json.dumps(payload["counts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
