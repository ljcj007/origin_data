#!/usr/bin/env python3
"""批量执行虚时统计、SAC 和谱平均，并合并为可绘图数据。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


def read_modes(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_beta(run_dir: Path) -> float:
    text = (run_dir / "sse_graph_summary.txt").read_text(encoding="utf-8")
    match = re.search(r"beta\s*=\s*([0-9Ee+\-.]+)", text)
    if not match:
        raise ValueError("无法从 sse_graph_summary.txt 读取 beta")
    return float(match.group(1))


def run_command(command: list[str], cwd: Path, name: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    (cwd / f"{name}.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (cwd / f"{name}.stderr.txt").write_text(result.stderr, encoding="utf-8")
    return result


def q_file(mode_id: int) -> str:
    if mode_id < 1 or mode_id > 999:
        raise ValueError("原始 tres 文件名仅支持 1..999")
    return f"q{mode_id:03d}.dat"


def normalized_endpoint(path: Path) -> float:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split()
    nt = int(header[1])
    values = [float(lines[index].split()[1]) for index in range(1, nt + 1)]
    return values[-1]


def path_distances(modes: list[dict[str, str]]) -> dict[int, float]:
    first_channel = modes[0]["channel"]
    q_rows = [row for row in modes if row["channel"] == first_channel]
    distance = 0.0
    result: dict[int, float] = {}
    previous: tuple[float, float] | None = None
    for row in q_rows:
        point = (float(row["qx"]), float(row["qy"]))
        if previous is not None:
            distance += math.dist(previous, point)
        result[int(row["q_index"])] = distance
        previous = point
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--tres-exe", type=Path, default=Path("build/sac_adapted/tres.exe"))
    parser.add_argument("--sac-exe", type=Path, default=Path("build/sac_adapted/sac.exe"))
    parser.add_argument("--avspec-exe", type=Path, default=Path("build/sac_adapted/avspec.exe"))
    parser.add_argument("--channels", default="total,edge,cls")
    parser.add_argument("--beta", type=float)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--rebin", type=int, default=1)
    parser.add_argument("--time-stride", type=int, default=1)
    parser.add_argument("--feature", choices=("free", "delta"), default="free")
    parser.add_argument("--delta-count", type=int, default=500)
    parser.add_argument("--output-spacing", type=float, default=0.02)
    parser.add_argument("--sampling-spacing", type=float, default=0.0005)
    parser.add_argument("--omega-max", type=float, default=4.0)
    parser.add_argument("--theta-steps", type=int, default=500)
    parser.add_argument("--spectrum-bins", type=int, default=50)
    parser.add_argument("--steps-per-bin", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=1675787505326236690)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"输出目录已存在：{output_dir}")
    output_dir.mkdir(parents=True)
    executables = {
        "tres": args.tres_exe.resolve(),
        "sac": args.sac_exe.resolve(),
        "avspec": args.avspec_exe.resolve(),
    }
    for name, path in executables.items():
        if not path.is_file():
            raise FileNotFoundError(f"缺少 {name}：{path}")

    modes = read_modes(run_dir / "dynamics_modes.csv")
    selected_channels = {value.strip() for value in args.channels.split(",") if value.strip()}
    selected = [row for row in modes if row["channel"] in selected_channels]
    if not selected:
        raise ValueError("没有匹配的动力学通道")
    beta = args.beta if args.beta is not None else read_beta(run_dir)
    nmode = len(modes)
    distances = path_distances(modes)
    records: list[dict[str, object]] = []
    spectra: list[dict[str, object]] = []

    for row in selected:
        mode_id = int(row["mode_id"])
        work = output_dir / f"mode-{mode_id:03d}-{row['channel']}"
        work.mkdir()
        shutil.copyfile(run_dir / "dynamics_cor.dat", work / "cor.dat")
        shutil.copyfile(run_dir / "dynamics_tgrid.dat", work / "tgrid.dat")
        (work / "tres.in").write_text(
            f"{nmode}\n{beta:.17g}\n{mode_id}\n0\n{args.rebin}\n{args.time_stride}\n{args.bootstrap}\n",
            encoding="utf-8",
        )
        (work / "seed.in").write_text(f"{args.seed + mode_id * 104729}\n", encoding="utf-8")
        tres_result = run_command([str(executables["tres"])], work, "tres")
        record: dict[str, object] = {
            "mode_id": mode_id,
            "q_index": int(row["q_index"]),
            "q_label": int(row["q_label"]),
            "channel": row["channel"],
            "tres_returncode": tres_result.returncode,
        }
        target = work / q_file(mode_id)
        if tres_result.returncode != 0 or not target.is_file():
            record["status"] = "tres_failed"
            records.append(record)
            continue
        endpoint = normalized_endpoint(target)
        record["normalized_endpoint"] = endpoint
        if not 0.0 < endpoint < 0.999:
            record["status"] = "nondecaying_or_invalid"
            records.append(record)
            continue

        shutil.copyfile(target, work / "t.in")
        feature = 0 if args.feature == "free" else 1
        (work / "sac.in").write_text(
            f"{feature} {args.delta_count} 1.0 {args.output_spacing:.17g}\n"
            f"{args.sampling_spacing:.17g} 0.0 {args.omega_max:.17g}\n"
            f"{args.theta_steps} {args.spectrum_bins} {args.steps_per_bin}\n",
            encoding="utf-8",
        )
        (work / "seed.in").write_text(f"{args.seed + mode_id * 130363}\n", encoding="utf-8")
        sac_result = run_command([str(executables["sac"])], work, "sac")
        record["sac_returncode"] = sac_result.returncode
        if sac_result.returncode != 0:
            record["status"] = "sac_failed"
            records.append(record)
            continue
        average_result = run_command([str(executables["avspec"])], work, "avspec")
        record["avspec_returncode"] = average_result.returncode
        spectrum_path = work / "sw.dat"
        if average_result.returncode != 0 or not spectrum_path.is_file():
            record["status"] = "average_failed"
            records.append(record)
            continue

        spectrum = np.loadtxt(spectrum_path)
        peak_index = int(np.argmax(spectrum[:, 1]))
        record["status"] = "complete"
        record["peak_omega"] = float(spectrum[peak_index, 0])
        record["spectral_integral"] = float(np.trapz(spectrum[:, 1], spectrum[:, 0]))
        records.append(record)
        for omega, sqw, detailed_balance in spectrum:
            spectra.append(
                {
                    "mode_id": mode_id,
                    "q_index": int(row["q_index"]),
                    "q_label": int(row["q_label"]),
                    "channel": row["channel"],
                    "qx": float(row["qx"]),
                    "qy": float(row["qy"]),
                    "path_distance": distances[int(row["q_index"])],
                    "omega": float(omega),
                    "S_q_omega": float(sqw),
                    "detailed_balance_spectrum": float(detailed_balance),
                }
            )

    summary = {
        "format": "sac-batch-summary-v1",
        "run_dir": str(run_dir),
        "beta": beta,
        "selected_channels": sorted(selected_channels),
        "settings": {
            "bootstrap": args.bootstrap,
            "rebin": args.rebin,
            "time_stride": args.time_stride,
            "feature": args.feature,
            "delta_count": args.delta_count,
            "output_spacing": args.output_spacing,
            "sampling_spacing": args.sampling_spacing,
            "omega_max": args.omega_max,
            "theta_steps": args.theta_steps,
            "spectrum_bins": args.spectrum_bins,
            "steps_per_bin": args.steps_per_bin,
        },
        "counts": {
            "selected": len(selected),
            "complete": sum(record["status"] == "complete" for record in records),
            "skipped_or_failed": sum(record["status"] != "complete" for record in records),
        },
        "records": records,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if spectra:
        with (output_dir / "spectra.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(spectra[0]))
            writer.writeheader()
            writer.writerows(spectra)
    print(json.dumps(summary["counts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
