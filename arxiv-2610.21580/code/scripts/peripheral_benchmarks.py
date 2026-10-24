#!/usr/bin/env python3
"""重建论文附录中的外围数值基准并记录来源边界。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "peripheral_benchmarks_v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def parse_summary(path: Path) -> dict[str, dict[str, float]]:
    values: dict[str, dict[str, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        tokens = line.split()
        if len(tokens) != 3:
            continue
        try:
            values[tokens[0]] = {"value": float(tokens[1]), "error": float(tokens[2])}
        except ValueError:
            continue
    return values


def polynomial_matrix_product(
    left: list[np.ndarray], right: list[np.ndarray], order: int
) -> list[np.ndarray]:
    result = [np.zeros_like(left[0]) for _ in range(order + 1)]
    for degree in range(order + 1):
        for first in range(degree + 1):
            result[degree] += left[first] @ right[degree - first]
    return result


def exact_square_ising_moments(size: int, temperature: float) -> dict[str, float]:
    """用传递矩阵的四阶 Taylor jet 求有限环面精确矩。"""
    if size < 2 or size % 2:
        raise ValueError("size 必须是至少为 2 的偶数")
    if temperature <= 0.0:
        raise ValueError("temperature 必须为正数")
    states = np.arange(1 << size, dtype=np.uint64)
    spins = np.empty((len(states), size), dtype=np.float64)
    for site in range(size):
        spins[:, site] = 2.0 * ((states >> np.uint64(site)) & 1) - 1.0
    row_energy = np.sum(spins * np.roll(spins, -1, axis=1), axis=1)
    row_magnetization = np.sum(spins, axis=1)
    vertical = spins @ spins.T
    beta = 1.0 / temperature
    exponent = beta * (
        vertical + 0.5 * (row_energy[:, None] + row_energy[None, :])
    )
    field_weight = 0.5 * (
        row_magnetization[:, None] + row_magnetization[None, :]
    )
    base = np.exp(exponent)
    transfer = [
        base * field_weight**degree / math.factorial(degree)
        for degree in range(5)
    ]
    identity = [np.zeros_like(base) for _ in range(5)]
    identity[0] = np.eye(len(states))
    power = identity
    factor = transfer
    exponent_size = size
    while exponent_size:
        if exponent_size & 1:
            power = polynomial_matrix_product(power, factor, 4)
        exponent_size //= 2
        if exponent_size:
            factor = polynomial_matrix_product(factor, factor, 4)
    coefficients = [float(np.trace(matrix)) for matrix in power]
    partition = coefficients[0]
    m2 = 2.0 * coefficients[2] / partition
    m4 = 24.0 * coefficients[4] / partition
    binder = 1.0 - m4 / (3.0 * m2 * m2)
    return {
        "size": size,
        "temperature": temperature,
        "partition_function": partition,
        "magnetization_m2": m2,
        "magnetization_m4": m4,
        "binder": binder,
    }


def compile_fortran(source: Path, executable: Path) -> dict[str, object]:
    executable.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "gfortran",
        "-std=f2008",
        "-Wall",
        "-Wextra",
        "-O2",
        "-J",
        str(executable.parent),
        str(source),
        "-o",
        str(executable),
    ]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"编译失败：{source}\n{completed.stderr}")
    return {
        "source": str(source.relative_to(ROOT).as_posix()),
        "source_sha256": sha256_file(source),
        "executable": str(executable.relative_to(ROOT).as_posix()),
        "executable_sha256": sha256_file(executable),
        "command": command,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def run_fortran_cell(
    executable: Path, cell_dir: Path, input_text: str, output_prefix: str
) -> dict[str, object]:
    cell_dir.mkdir(parents=True)
    input_path = cell_dir / "run.in"
    input_path.write_text(input_text, encoding="utf-8")
    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable.resolve()), str(input_path.resolve())],
        cwd=cell_dir,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    elapsed = time.perf_counter() - started
    (cell_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (cell_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    bins = cell_dir / f"{output_prefix}_bins.csv"
    summary = cell_dir / f"{output_prefix}_summary.txt"
    if completed.returncode != 0 or not bins.is_file() or not summary.is_file():
        raise RuntimeError(f"运行失败：{cell_dir}\n{completed.stderr}")
    return {
        "returncode": completed.returncode,
        "elapsed_seconds": elapsed,
        "input_sha256": sha256_file(input_path),
        "bins_sha256": sha256_file(bins),
        "summary_sha256": sha256_file(summary),
        "observables": parse_summary(summary),
    }


def exact_heisenberg_square(graph: Path) -> dict[str, object]:
    from ed_graph_heisenberg import lowest_state, make_lattice

    nsite, bonds, _, _, _ = make_lattice(graph, 4, 4)
    energy, _, states, residual, norm = lowest_state(
        nsite, bonds, nsite // 2, 1.0e-12, 1000, False
    )
    return {
        "method": "scipy.sparse.linalg.eigsh in Sz=0",
        "nsite": nsite,
        "nbond": len(bonds),
        "sector_dimension": len(states),
        "ground_state_energy": energy,
        "ground_state_energy_per_site": energy / nsite,
        "residual_norm": residual,
        "state_norm": norm,
    }


def exact_heisenberg_square_octagon() -> dict[str, object]:
    from ed_square_octagon import lowest_state, make_hamiltonian, make_lattice

    nsite, bonds, _ = make_lattice(2, 2)
    matrix, _ = make_hamiltonian(nsite, bonds, nsite // 2)
    energy, _, residual = lowest_state(matrix, maxiter=160)
    return {
        "method": "fully reorthogonalized Lanczos in Sz=0",
        "nsite": nsite,
        "nbond": len(bonds),
        "sector_dimension": len(matrix[0]),
        "ground_state_energy": energy,
        "ground_state_energy_per_site": energy / nsite,
        "reported_lanczos_residual": residual,
    }


def standardized_difference(value: float, error: float, oracle: float, oracle_error: float) -> float:
    denominator = math.hypot(error, oracle_error)
    if denominator == 0.0:
        return math.inf if value != oracle else 0.0
    return (value - oracle) / denominator


def build_comparisons(
    runs: dict[str, dict[str, object]],
    ising_exact: dict[str, float],
    square_ed: dict[str, object],
    octagon_ed: dict[str, object],
) -> list[dict[str, object]]:
    declared = {
        "square_ising_l8": (0.6137, 0.0026),
        "square_heisenberg_l4": (-0.701696, 0.000261),
        "square_heisenberg_l8": (-0.673329, 0.000211),
        "square_octagon_n16": (-0.570028, 0.000171),
    }
    specifications = [
        (
            "square_ising_l8",
            "binder_staggered",
            float(ising_exact["binder"]),
            0.0,
            "finite L=8 transfer matrix",
            "local_exact",
        ),
        (
            "square_heisenberg_l4",
            "energy_per_site",
            float(square_ed["ground_state_energy_per_site"]),
            float(square_ed["residual_norm"]) / int(square_ed["nsite"]),
            "N=16 sparse exact diagonalization",
            "local_exact",
        ),
        (
            "square_heisenberg_l8",
            "energy_per_site",
            -0.673487,
            0.000004,
            "Sandvik 1997 value quoted by manuscript",
            "primary_source_not_available_offline",
        ),
        (
            "square_octagon_n16",
            "energy_per_site",
            float(octagon_ed["ground_state_energy_per_site"]),
            0.0,
            "N=16 fully reorthogonalized Lanczos",
            "local_exact",
        ),
    ]
    rows: list[dict[str, object]] = []
    for name, observable, oracle, oracle_error, oracle_method, source_status in specifications:
        measured = runs[name]["observables"][observable]
        value = float(measured["value"])
        error = float(measured["error"])
        quoted_value, quoted_error = declared[name]
        rows.append(
            {
                "benchmark": name,
                "observable": observable,
                "rerun_value": value,
                "rerun_standard_error": error,
                "manuscript_value": quoted_value,
                "manuscript_standard_error": quoted_error,
                "rerun_vs_manuscript_z": standardized_difference(
                    value, error, quoted_value, quoted_error
                ),
                "oracle_value": oracle,
                "oracle_uncertainty": oracle_error,
                "oracle_method": oracle_method,
                "oracle_source_status": source_status,
                "rerun_vs_oracle_z": standardized_difference(
                    value, error, oracle, oracle_error
                ),
                "cross_method_tag": "not_assigned_missing_user_accuracy_budget",
            }
        )
    return rows


def write_comparison_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_stiffness_checks(runs: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    specifications = {
        "square_heisenberg_l4": (4, 32.0),
        "square_heisenberg_l8": (8, 64.0),
    }
    rows = []
    for name, (size, beta) in specifications.items():
        observables = runs[name]["observables"]
        winding = float(observables["winding_squared"]["value"])
        stiffness = float(observables["spin_stiffness_cell_average"]["value"])
        expected = 3.0 * winding / (4.0 * beta)
        rows.append(
            {
                "benchmark": name,
                "L": size,
                "beta": beta,
                "winding_squared": winding,
                "stiffness": stiffness,
                "stiffness_from_winding": expected,
                "identity_residual": stiffness - expected,
                "stiffness_x": float(observables["spin_stiffness_cell_x"]["value"]),
                "stiffness_y": float(observables["spin_stiffness_cell_y"]["value"]),
            }
        )
    return rows


def copy_legacy_snapshots(output: Path) -> None:
    legacy_dir = output / "legacy"
    legacy_dir.mkdir(exist_ok=True)
    groups = {
        "square_ising": ["ising_graph.in", "ising_bins.csv", "ising_summary.txt"],
        "square_octagon": [
            "sse_graph.in",
            "sse_graph_bins.csv",
            "sse_graph_summary.txt",
        ],
    }
    records: dict[str, object] = {
        "status": "legacy_unmanifested",
        "boundary": (
            "These root-level artifacts predate formal manifests. Their outputs are preserved, "
            "but no contemporaneous source/executable hash was recorded."
        ),
        "groups": {},
    }
    for group, names in groups.items():
        destination = legacy_dir / group
        destination.mkdir(exist_ok=True)
        files = []
        for name in names:
            source = ROOT / name
            target = destination / name
            shutil.copy2(source, target)
            files.append(
                {
                    "path": str(target.relative_to(output).as_posix()),
                    "sha256": sha256_file(target),
                    "size": target.stat().st_size,
                }
            )
        records["groups"][group] = files
    write_json(legacy_dir / "provenance_boundary.json", records)


def run(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"输出已存在，拒绝覆盖：{output}")
    output.mkdir(parents=True)
    graph_dir = output / "graphs"
    graph_dir.mkdir()
    square_graph = graph_dir / "square.dat"
    octagon_graph = graph_dir / "square_octagon.dat"
    shutil.copy2(ROOT / "graphs" / "square.dat", square_graph)
    shutil.copy2(ROOT / "graphs" / "square_octagon_manual.dat", octagon_graph)

    compiler = subprocess.run(
        ["gfortran", "--version"], text=True, capture_output=True, check=True
    ).stdout.splitlines()[0]
    ising_executable = output / "build" / "ising_graph.exe"
    sse_executable = output / "build" / "sse_graph.exe"
    provenance = {
        "python": sys.version,
        "numpy": np.__version__,
        "compiler": compiler,
        "ising": compile_fortran(ROOT / "ising_graph.f90", ising_executable),
        "sse": compile_fortran(ROOT / "sse_graph.f90", sse_executable),
        "graphs": {
            "square_sha256": sha256_file(square_graph),
            "square_octagon_sha256": sha256_file(octagon_graph),
        },
        "analysis_sources": {
            "peripheral_benchmarks.py": sha256_file(Path(__file__)),
            "ed_graph_heisenberg.py": sha256_file(ROOT / "scripts" / "ed_graph_heisenberg.py"),
            "ed_square_octagon.py": sha256_file(ROOT / "scripts" / "ed_square_octagon.py"),
        },
    }
    write_json(output / "provenance.json", provenance)

    runs: dict[str, dict[str, object]] = {}
    runs["square_ising_l8"] = run_fortran_cell(
        ising_executable,
        output / "cells" / "square_ising_l8",
        f"{square_graph.resolve()}\n8 8\n2.269185314213022 0 1\n"
        "12 4000 2000 20260903\nmetropolis\n",
        "ising",
    )
    runs["square_heisenberg_l4"] = run_fortran_cell(
        sse_executable,
        output / "cells" / "square_heisenberg_l4",
        f"{square_graph.resolve()}\n4 4\n32 1\n1\n"
        "16 5000 5000 20260905\n",
        "sse_graph",
    )
    runs["square_heisenberg_l8"] = run_fortran_cell(
        sse_executable,
        output / "cells" / "square_heisenberg_l8",
        f"{square_graph.resolve()}\n8 8\n64 1\n1\n"
        "12 3000 3000 20260906\n",
        "sse_graph",
    )
    runs["square_octagon_n16"] = run_fortran_cell(
        sse_executable,
        output / "cells" / "square_octagon_n16",
        f"{octagon_graph.resolve()}\n2 2\n16 2\n1 1\n"
        "16 10000 10000 20260903\n",
        "sse_graph",
    )
    for name, payload in runs.items():
        write_json(output / "cells" / name / "manifest.json", payload)

    exact_dir = output / "exact"
    exact_dir.mkdir()
    ising_exact = exact_square_ising_moments(8, 2.269185314213022)
    square_ed = exact_heisenberg_square(square_graph)
    octagon_ed = exact_heisenberg_square_octagon()
    write_json(exact_dir / "square_ising_l8_transfer.json", ising_exact)
    write_json(exact_dir / "square_heisenberg_l4_ed.json", square_ed)
    write_json(exact_dir / "square_octagon_n16_ed.json", octagon_ed)

    rows = build_comparisons(runs, ising_exact, square_ed, octagon_ed)
    stiffness_checks = build_stiffness_checks(runs)
    write_comparison_csv(output / "cross-method-check.csv", rows)
    local_exact_rows = [row for row in rows if row["oracle_source_status"] == "local_exact"]
    summary = {
        "schema_version": 1,
        "status": "partial_offline_closure",
        "benchmarks_reconstructed": len(rows),
        "benchmarks_with_local_independent_oracle": len(local_exact_rows),
        "maximum_abs_z_against_local_oracle": max(
            abs(float(row["rerun_vs_oracle_z"])) for row in local_exact_rows
        ),
        "maximum_abs_z_against_manuscript_numbers": max(
            abs(float(row["rerun_vs_manuscript_z"])) for row in rows
        ),
        "blocking_gap": (
            "L=8 square-Heisenberg reference -0.673487(4) is only present as a "
            "manuscript-declared Sandvik-1997 value; the primary source is absent offline."
        ),
        "accuracy_budget": (
            "No Agreement/Disagreement tag is assigned because the user has not supplied "
            "cross-method accuracy budgets; signed combined-standard-error distances are recorded."
        ),
        "stiffness_normalization": {
            "formula": "rho_s = 3 <W_x^2 + W_y^2> / (4 beta)",
            "status": "nonzero_internal_identity_pass; external thermodynamic benchmark pending",
            "maximum_absolute_identity_residual": max(
                abs(float(row["identity_residual"])) for row in stiffness_checks
            ),
            "checks": stiffness_checks,
        },
        "comparisons": rows,
    }
    write_json(output / "validation_summary.json", summary)
    copy_legacy_snapshots(output)


def resume_analysis(output: Path) -> None:
    if not output.is_dir():
        raise FileNotFoundError(f"结果目录不存在：{output}")
    required = [
        "square_ising_l8",
        "square_heisenberg_l4",
        "square_heisenberg_l8",
        "square_octagon_n16",
    ]
    runs = {
        name: json.loads(
            (output / "cells" / name / "manifest.json").read_text(encoding="utf-8")
        )
        for name in required
    }
    provenance_path = output / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["analysis_python"] = sys.version
    provenance["analysis_numpy"] = np.__version__
    provenance["analysis_sources"]["peripheral_benchmarks.py"] = sha256_file(Path(__file__))
    write_json(provenance_path, provenance)

    exact_dir = output / "exact"
    exact_dir.mkdir(exist_ok=True)
    ising_exact = exact_square_ising_moments(8, 2.269185314213022)
    square_ed = exact_heisenberg_square(output / "graphs" / "square.dat")
    octagon_ed = exact_heisenberg_square_octagon()
    write_json(exact_dir / "square_ising_l8_transfer.json", ising_exact)
    write_json(exact_dir / "square_heisenberg_l4_ed.json", square_ed)
    write_json(exact_dir / "square_octagon_n16_ed.json", octagon_ed)

    rows = build_comparisons(runs, ising_exact, square_ed, octagon_ed)
    stiffness_checks = build_stiffness_checks(runs)
    write_comparison_csv(output / "cross-method-check.csv", rows)
    local_exact_rows = [row for row in rows if row["oracle_source_status"] == "local_exact"]
    summary = {
        "schema_version": 1,
        "status": "partial_offline_closure",
        "benchmarks_reconstructed": len(rows),
        "benchmarks_with_local_independent_oracle": len(local_exact_rows),
        "maximum_abs_z_against_local_oracle": max(
            abs(float(row["rerun_vs_oracle_z"])) for row in local_exact_rows
        ),
        "maximum_abs_z_against_manuscript_numbers": max(
            abs(float(row["rerun_vs_manuscript_z"])) for row in rows
        ),
        "blocking_gap": (
            "L=8 square-Heisenberg reference -0.673487(4) is only present as a "
            "manuscript-declared Sandvik-1997 value; the primary source is absent offline."
        ),
        "accuracy_budget": (
            "No Agreement/Disagreement tag is assigned because the user has not supplied "
            "cross-method accuracy budgets; signed combined-standard-error distances are recorded."
        ),
        "stiffness_normalization": {
            "formula": "rho_s = 3 <W_x^2 + W_y^2> / (4 beta)",
            "status": "nonzero_internal_identity_pass; external thermodynamic benchmark pending",
            "maximum_absolute_identity_residual": max(
                abs(float(row["identity_residual"])) for row in stiffness_checks
            ),
            "checks": stiffness_checks,
        },
        "comparisons": rows,
    }
    write_json(output / "validation_summary.json", summary)
    copy_legacy_snapshots(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume-analysis", action="store_true")
    args = parser.parse_args()
    if args.resume_analysis:
        resume_analysis(args.output.resolve())
    else:
        run(args.output.resolve())


if __name__ == "__main__":
    main()
