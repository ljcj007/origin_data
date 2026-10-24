#!/usr/bin/env python3
"""用独立二聚体精确解校验 SSE 输出。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


LINEAR_GATE = (
    "energy_per_site",
    "staggered_m2",
    "staggered_m4",
    "ferrimagnetic_m2",
    "sublattice_plus_m2",
    "sublattice_minus_m2",
    "sublattice_cross",
    "mean_expansion_order",
)
ZSCORE_LIMIT = 4.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def mean_stderr(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        raise ValueError("至少需要两个 bin")
    mean = math.fsum(values) / len(values)
    error = math.sqrt(
        math.fsum((value - mean) ** 2 for value in values)
        / (len(values) * (len(values) - 1))
    )
    return mean, error


def blocked_stderr(values: list[float], block_size: int) -> float:
    if len(values) % block_size != 0:
        raise ValueError("bin 数不能被重分块长度整除")
    blocked = [
        math.fsum(values[start : start + block_size]) / block_size
        for start in range(0, len(values), block_size)
    ]
    return mean_stderr(blocked)[1]


def lag_one(values: list[float]) -> float:
    mean = math.fsum(values) / len(values)
    denominator = math.fsum((value - mean) ** 2 for value in values)
    if denominator == 0.0:
        return 0.0
    return math.fsum(
        (values[index] - mean) * (values[index + 1] - mean)
        for index in range(len(values) - 1)
    ) / denominator


def exact_dimer_observables(beta: float, coupling: float, ndimer: int) -> dict[str, float]:
    """返回 ndimer 个独立自旋 1/2 二聚体的程序归一化精确值。"""
    x = math.exp(beta * coupling)
    denominator = x + 3.0
    nsite = 2 * ndimer
    energy_dimer = 0.75 * coupling * (1.0 - x) / denominator
    staggered_single_second = (x + 1.0) / denominator
    staggered_total_second = ndimer * staggered_single_second
    staggered_total_fourth = (
        ndimer * staggered_single_second
        + 3.0 * ndimer * (ndimer - 1) * staggered_single_second**2
    )
    total_spin_squared = 6.0 * ndimer / denominator
    mean_n = beta * coupling * ndimer * x / denominator
    return {
        "energy_per_site": energy_dimer / 2.0,
        "specific_heat_per_site": 1.5 * (beta * coupling) ** 2 * x / denominator**2,
        "staggered_m2": 3.0 * staggered_total_second / nsite**2,
        "staggered_m4": staggered_total_fourth / nsite**4,
        "binder_z_staggered": 1.0
        - staggered_total_fourth / (3.0 * staggered_total_second**2),
        "ferrimagnetic_m2": total_spin_squared / nsite**2,
        "uniform_susceptibility": beta / denominator,
        "sublattice_plus_m2": 3.0 / (4.0 * ndimer),
        "sublattice_minus_m2": 3.0 / (4.0 * ndimer),
        "sublattice_cross": energy_dimer / (coupling * ndimer),
        "spin_stiffness_cell_x": 0.0,
        "spin_stiffness_cell_y": 0.0,
        "spin_stiffness_cell_average": 0.0,
        "spin_stiffness_times_L": 0.0,
        "winding_squared": 0.0,
        "mean_expansion_order": mean_n,
        "total_spin_squared": total_spin_squared,
    }


def read_input(run_dir: Path) -> dict[str, object]:
    lines = [line.strip() for line in (run_dir / "sse_graph.in").read_text().splitlines() if line.strip()]
    if len(lines) != 5:
        raise ValueError("SSE 输入应恰有五个非空行")
    lx, ly = (int(value) for value in lines[1].split())
    beta_text, ntypes_text = lines[2].split()
    coupling_values = [float(value) for value in lines[3].split()]
    nbins, equil_sweeps, bin_sweeps, seed = (int(value) for value in lines[4].split())
    ntypes = int(ntypes_text)
    if ntypes != 1 or len(coupling_values) != 1:
        raise ValueError("二聚体校验只接受一种交换常数")
    return {
        "graph": lines[0],
        "lx": lx,
        "ly": ly,
        "beta": float(beta_text),
        "ntypes": ntypes,
        "coupling": coupling_values[0],
        "nbins": nbins,
        "equil_sweeps": equil_sweeps,
        "bin_sweeps": bin_sweeps,
        "seed": seed,
    }


def validate_graph(path: Path) -> None:
    tokens = path.read_text().split()
    expected = ["2", "1", "3", "1", "2", "0", "0", "0", "1", "1", "-1", "0", "0", "0"]
    if tokens != expected:
        raise ValueError("图不是一条胞内键连接的二分二聚体")


def read_bins(path: Path, nbins: int) -> tuple[list[str], dict[str, list[float]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != nbins:
        raise ValueError(f"bin 数错误：{len(rows)} != {nbins}")
    if [int(row["bin"]) for row in rows] != list(range(1, nbins + 1)):
        raise ValueError("bin 编号不连续")
    names = [name for name in rows[0] if name != "bin"]
    columns = {name: [float(row[name]) for row in rows] for name in names}
    if not all(math.isfinite(value) for values in columns.values() for value in values):
        raise ValueError("原始 bin 含非有限值")
    return names, columns


def read_summary(path: Path) -> dict[str, tuple[float, float]]:
    result: dict[str, tuple[float, float]] = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            result[match.group(1)] = (float(match.group(2)), float(match.group(3)))
    return result


def zscore(mean: float, exact: float, stderr: float) -> float | None:
    if stderr == 0.0:
        return 0.0 if abs(mean - exact) <= 1.0e-14 else None
    return (mean - exact) / stderr


def max_identity_residual(columns: dict[str, list[float]], beta: float, nsite: int, ndimer: int) -> dict[str, float]:
    residuals = {"energy_expansion": 0.0, "susceptibility_ferrimagnetic": 0.0, "staggered_sublattice": 0.0}
    for index in range(len(columns["energy_per_site"])):
        energy_from_n = -columns["mean_expansion_order"][index] / (beta * nsite) + ndimer / (4.0 * nsite)
        residuals["energy_expansion"] = max(residuals["energy_expansion"], abs(columns["energy_per_site"][index] - energy_from_n))
        chi_from_m2 = beta * nsite * columns["ferrimagnetic_m2"][index] / 3.0
        residuals["susceptibility_ferrimagnetic"] = max(
            residuals["susceptibility_ferrimagnetic"],
            abs(columns["uniform_susceptibility"][index] - chi_from_m2),
        )
        staggered_from_sub = 3.0 * (
            columns["sublattice_plus_m2"][index] * ndimer**2 / 3.0
            + columns["sublattice_minus_m2"][index] * ndimer**2 / 3.0
            - 2.0 * columns["sublattice_cross"][index] * ndimer**2 / 3.0
        ) / nsite**2
        residuals["staggered_sublattice"] = max(
            residuals["staggered_sublattice"],
            abs(columns["staggered_m2"][index] - staggered_from_sub),
        )
    return residuals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    config = read_input(run_dir)
    graph_path = (run_dir / str(config["graph"])).resolve()
    validate_graph(graph_path)
    lx, ly = int(config["lx"]), int(config["ly"])
    ndimer = lx * ly
    nsite = 2 * ndimer
    exact = exact_dimer_observables(float(config["beta"]), float(config["coupling"]), ndimer)
    names, columns = read_bins(run_dir / "sse_graph_bins.csv", int(config["nbins"]))
    summary = read_summary(run_dir / "sse_graph_summary.txt")

    comparisons: dict[str, dict[str, object]] = {}
    for name in names:
        mean, naive_stderr = mean_stderr(columns[name])
        block_errors = {
            str(block_size): blocked_stderr(columns[name], block_size)
            for block_size in (1, 2, 4, 8)
        }
        stderr = max(block_errors.values())
        exact_value = exact[name]
        z = zscore(mean, exact_value, stderr)
        comparisons[name] = {
            "mean": mean,
            "naive_stderr": naive_stderr,
            "conservative_stderr": stderr,
            "blocked_stderr": block_errors,
            "lag_one_bin": lag_one(columns[name]),
            "exact": exact_value,
            "difference": mean - exact_value,
            "z_score": z,
            "statistical_gate": name in LINEAR_GATE,
        }
        if name not in summary:
            raise ValueError(f"summary 缺少 {name}")
        if abs(summary[name][0] - mean) > 5.0e-14 or abs(summary[name][1] - naive_stderr) > 5.0e-14:
            raise ValueError(f"summary 与 bins 不一致：{name}")

    stot_mean = comparisons["ferrimagnetic_m2"]["mean"] * nsite**2
    stot_stderr = comparisons["ferrimagnetic_m2"]["conservative_stderr"] * nsite**2
    stot_z = zscore(float(stot_mean), exact["total_spin_squared"], float(stot_stderr))
    gate_zscores = {name: comparisons[name]["z_score"] for name in LINEAR_GATE}
    stochastic_pass = all(value is not None and abs(float(value)) <= ZSCORE_LIMIT for value in gate_zscores.values())
    invariant_residuals = max_identity_residual(columns, float(config["beta"]), nsite, ndimer)
    invariant_pass = max(invariant_residuals.values()) <= 5.0e-14
    zero_shift_pass = all(
        all(abs(value) <= 1.0e-14 for value in columns[name])
        for name in (
            "spin_stiffness_cell_x",
            "spin_stiffness_cell_y",
            "spin_stiffness_cell_average",
            "spin_stiffness_times_L",
            "winding_squared",
        )
    )
    payload = {
        "schema_version": 1,
        "verdict": "pass" if stochastic_pass and invariant_pass and zero_shift_pass else "fail",
        "model": {
            "hamiltonian": "H = J * sum_d S_(d,1) dot S_(d,2)",
            "spin": 0.5,
            "coupling": config["coupling"],
            "beta": config["beta"],
            "lx": lx,
            "ly": ly,
            "independent_dimers": ndimer,
            "nsite": nsite,
            "levels_per_dimer": {"singlet": -0.75 * float(config["coupling"]), "triplet": 0.25 * float(config["coupling"]), "triplet_degeneracy": 3},
        },
        "sampling": {
            "nbins": config["nbins"],
            "equil_sweeps": config["equil_sweeps"],
            "bin_sweeps": config["bin_sweeps"],
            "seed": config["seed"],
        },
        "statistical_rule": {"gate_observables": list(LINEAR_GATE), "maximum_absolute_z_score": ZSCORE_LIMIT},
        "comparisons": comparisons,
        "total_spin_squared": {"mean": stot_mean, "stderr": stot_stderr, "exact": exact["total_spin_squared"], "z_score": stot_z},
        "identity_max_absolute_residual": invariant_residuals,
        "zero_shift_stiffness_exactly_zero": zero_shift_pass,
        "provenance": {
            "validator": Path(__file__).resolve().as_posix(),
            "validator_sha256": sha256(Path(__file__).resolve()),
            "source": args.source.resolve().as_posix(),
            "source_sha256": sha256(args.source.resolve()),
            "executable": args.executable.resolve().as_posix(),
            "executable_sha256": sha256(args.executable.resolve()),
            "graph_sha256": sha256(graph_path),
            "input_sha256": sha256(run_dir / "sse_graph.in"),
        },
    }
    output = (args.output or run_dir / "exact_comparison.json").resolve()
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"二聚体 SSE 校验：{payload['verdict']}；最大 |z|={max(abs(float(value)) for value in gate_zscores.values() if value is not None):.3f}")
    if payload["verdict"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
