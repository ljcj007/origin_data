#!/usr/bin/env python3
"""汇总 CLS 谱的尺寸与解析延拓稳定性。"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--case", nargs=3, action="append", metavar=("LABEL", "L", "DIRECTORY"), required=True
    )
    args = parser.parse_args()

    records: list[dict[str, object]] = []
    cases = []
    for label, size_text, directory_text in args.case:
        size = int(size_text)
        directory = Path(directory_text).resolve()
        rows = [row for row in read_rows(directory / "spectra.csv") if row["channel"] == "cls"]
        q_indices = sorted({int(row["q_index"]) for row in rows})
        case_records = []
        for q_index in q_indices:
            selected = [row for row in rows if int(row["q_index"]) == q_index]
            omega = np.asarray([float(row["omega"]) for row in selected])
            intensity = np.asarray([float(row["S_q_omega"]) for row in selected])
            integral = float(np.trapz(intensity, omega))
            centroid = float(np.trapz(omega * intensity, omega) / integral)
            variance = float(np.trapz((omega - centroid) ** 2 * intensity, omega) / integral)
            record = {
                "case": label,
                "L": size,
                "q_index": q_index,
                "path_distance": float(selected[0]["path_distance"]),
                "peak_omega": float(omega[int(np.argmax(intensity))]),
                "centroid_omega": centroid,
                "rms_width": variance**0.5,
                "spectral_integral": integral,
            }
            case_records.append(record)
            records.append(record)
        peaks = [float(record["peak_omega"]) for record in case_records]
        centroids = [float(record["centroid_omega"]) for record in case_records]
        integrals = [float(record["spectral_integral"]) for record in case_records]
        cases.append(
            {
                "case": label,
                "L": size,
                "directory": str(directory),
                "q_count": len(case_records),
                "peak_mean": float(np.mean(peaks)),
                "peak_min": min(peaks),
                "peak_max": max(peaks),
                "peak_bandwidth": max(peaks) - min(peaks),
                "centroid_mean": float(np.mean(centroids)),
                "centroid_min": min(centroids),
                "centroid_max": max(centroids),
                "centroid_bandwidth": max(centroids) - min(centroids),
                "integral_mean": float(np.mean(integrals)),
                "integral_relative_range": (max(integrals) - min(integrals)) / float(np.mean(integrals)),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.with_suffix(".csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    payload = {"format": "sac-stability-v1", "cases": cases}
    args.output.with_suffix(".json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(cases, ensure_ascii=False))


if __name__ == "__main__":
    main()
