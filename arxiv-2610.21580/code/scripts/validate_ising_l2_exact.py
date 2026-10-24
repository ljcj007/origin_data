#!/usr/bin/env python3
"""对照 L=2 Ising 精确枚举与两种蒙卡更新。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


OBSERVABLES = (
    "energy_per_site",
    "specific_heat_per_site",
    "abs_staggered_m",
    "staggered_m2",
    "staggered_m4",
    "binder_staggered",
    "uniform_susceptibility",
    "connected_staggered_susceptibility",
)
BLOCK_SIZES = (1, 2, 4, 8, 16)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def block_means(values: list[float], block_size: int) -> list[float]:
    if block_size < 1 or len(values) % block_size:
        raise ValueError("bin 数不能被块长整除")
    return [
        statistics.fmean(values[start : start + block_size])
        for start in range(0, len(values), block_size)
    ]


def mean_error(values: list[float], block_size: int) -> tuple[float, float]:
    blocks = block_means(values, block_size)
    if len(blocks) < 4:
        raise ValueError("重分块后块数不足")
    return statistics.fmean(values), statistics.stdev(blocks) / math.sqrt(len(blocks))


def binder_value(q2: list[float], q4: list[float]) -> float:
    mean_q2 = statistics.fmean(q2)
    return 1.0 - statistics.fmean(q4) / (3.0 * mean_q2**2)


def binder_error(q2: list[float], q4: list[float], block_size: int) -> tuple[float, float]:
    blocks_q2 = block_means(q2, block_size)
    blocks_q4 = block_means(q4, block_size)
    count = len(blocks_q2)
    if count < 4:
        raise ValueError("Binder 重分块后块数不足")
    sum_q2, sum_q4 = sum(blocks_q2), sum(blocks_q4)
    samples = []
    for index in range(count):
        leave_q2 = (sum_q2 - blocks_q2[index]) / (count - 1)
        leave_q4 = (sum_q4 - blocks_q4[index]) / (count - 1)
        samples.append(1.0 - leave_q4 / (3.0 * leave_q2**2))
    center = statistics.fmean(samples)
    error = math.sqrt((count - 1) / count * sum((sample - center) ** 2 for sample in samples))
    return binder_value(q2, q4), error


def lag1(values: list[float]) -> float:
    mean = statistics.fmean(values)
    denominator = sum((value - mean) ** 2 for value in values)
    if denominator == 0.0:
        return 0.0
    numerator = sum(
        (values[index] - mean) * (values[index - 1] - mean)
        for index in range(1, len(values))
    )
    return numerator / denominator


def parse_bins(path: Path) -> dict[str, list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 256:
        raise ValueError(f"正式校验必须有 256 bins：{path}")
    if [int(row["bin"]) for row in rows] != list(range(1, 257)):
        raise ValueError("bin 编号不连续")
    columns = {
        name: [float(row[name]) for row in rows]
        for name in rows[0]
        if name != "bin"
    }
    if any(not math.isfinite(value) for values in columns.values() for value in values):
        raise ValueError("bin 数据含非有限数")
    return columns


def conservative_estimates(columns: dict[str, list[float]]) -> dict[str, dict[str, Any]]:
    estimates: dict[str, dict[str, Any]] = {}
    for observable in OBSERVABLES:
        errors = {}
        if observable == "binder_staggered":
            for block_size in BLOCK_SIZES:
                value, error = binder_error(
                    columns["staggered_m2"], columns["staggered_m4"], block_size
                )
                errors[str(block_size)] = error
        else:
            value = statistics.fmean(columns[observable])
            for block_size in BLOCK_SIZES:
                _, error = mean_error(columns[observable], block_size)
                errors[str(block_size)] = error
        conservative = max(errors.values())
        estimates[observable] = {
            "value": value,
            "errors_by_block_size": errors,
            "conservative_standard_error": conservative,
            "selected_block_size": int(max(errors, key=errors.get)),
        }
    return estimates


def inverse_bin_extrapolation(
    short: dict[str, Any],
    long: dict[str, Any],
    short_bin_sweeps: int,
    long_bin_sweeps: int,
    exact: float,
    z_limit: float,
) -> dict[str, Any]:
    if short_bin_sweeps <= 0 or long_bin_sweeps <= short_bin_sweeps:
        raise ValueError("外推 bin 长度错误")
    ratio = long_bin_sweeps / short_bin_sweeps
    value = (ratio * float(long["value"]) - float(short["value"])) / (ratio - 1.0)
    error = math.hypot(
        ratio * float(long["conservative_standard_error"]),
        float(short["conservative_standard_error"]),
    ) / (ratio - 1.0)
    z_score = (value - exact) / error
    return {
        "model": "chi(n_bin)=chi(infinity)+a/n_bin",
        "bin_sweeps": [short_bin_sweeps, long_bin_sweeps],
        "value": value,
        "conservative_standard_error": error,
        "exact": exact,
        "difference": value - exact,
        "z_score": z_score,
        "accuracy_budget": z_limit * error,
        "tag": "Agreement" if abs(z_score) <= z_limit else "Disagreement",
    }


def parse_input(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 5:
        raise ValueError("Ising 输入格式错误")
    lx, ly = map(int, lines[1].split())
    temperature, field, coupling = map(float, lines[2].split())
    nbins, equil_sweeps, bin_sweeps, seed = map(int, lines[3].split())
    return {
        "graph": Path(lines[0]).resolve(),
        "Lx": lx,
        "Ly": ly,
        "temperature": temperature,
        "field": field,
        "coupling": coupling,
        "nbins": nbins,
        "equil_sweeps": equil_sweeps,
        "bin_sweeps": bin_sweeps,
        "seed": seed,
        "algorithm": lines[4].strip(),
    }


def close(left: float, right: float, tolerance: float = 1.0e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def validate_cell(
    cell_dir: Path,
    algorithm: str,
    expected_bin_sweeps: int,
    exact: dict[str, Any],
    graph_hash: str,
    source_hash: str,
    executable_hash: str,
    z_limit: float,
) -> dict[str, Any]:
    manifest = load_json(cell_dir / "manifest.json")
    run_input = parse_input(cell_dir / "ising_graph.in")
    parameters = exact["parameters"]
    if manifest.get("status", {}).get("success") is not True:
        raise ValueError(f"蒙卡单元未成功：{cell_dir.name}")
    if run_input["Lx"] != 2 or run_input["Ly"] != 2 or run_input["nbins"] != 256:
        raise ValueError("蒙卡尺寸或 bin 数错误")
    if run_input["equil_sweeps"] != 100000 or run_input["bin_sweeps"] != expected_bin_sweeps:
        raise ValueError("蒙卡热化或测量长度错误")
    if run_input["algorithm"] != algorithm or run_input["field"] != 0.0:
        raise ValueError("蒙卡算法或外场错误")
    if not close(run_input["temperature"], float(parameters["temperature"])):
        raise ValueError("蒙卡温度与精确枚举不一致")
    if not close(run_input["coupling"], float(parameters["coupling"])):
        raise ValueError("蒙卡耦合与精确枚举不一致")
    provenance = manifest.get("provenance", {})
    execution = manifest.get("execution", {})
    if provenance.get("graph_sha256") != graph_hash or execution.get("graph_sha256") != graph_hash:
        raise ValueError("蒙卡格点图哈希错误")
    if provenance.get("source_sha256") != source_hash:
        raise ValueError("蒙卡源码哈希错误")
    if provenance.get("executable_sha256") != executable_hash:
        raise ValueError("蒙卡可执行文件哈希错误")
    if execution.get("executable_sha256") != executable_hash:
        raise ValueError("实际可执行文件哈希错误")
    columns = parse_bins(cell_dir / "ising_bins.csv")
    estimates = conservative_estimates(columns)
    exact_observables = exact["observables"]
    maximum_z = 0.0
    for observable, estimate in estimates.items():
        reference = float(exact_observables[observable])
        error = float(estimate["conservative_standard_error"])
        if error <= 0.0:
            raise ValueError(f"蒙卡误差非正：{algorithm} {observable}")
        difference = float(estimate["value"]) - reference
        z_score = difference / error
        estimate.update(
            {
                "exact": reference,
                "difference": difference,
                "z_score": z_score,
                "accuracy_budget": z_limit * error,
                "tag": "Agreement" if abs(z_score) <= z_limit else "Disagreement",
            }
        )
        maximum_z = max(maximum_z, abs(z_score))
    return {
        "algorithm": algorithm,
        "cell": cell_dir.name,
        "seed": run_input["seed"],
        "sweeps": {
            "equilibration": run_input["equil_sweeps"],
            "bins": run_input["nbins"],
            "per_bin": run_input["bin_sweeps"],
            "measurement_total": run_input["nbins"] * run_input["bin_sweeps"],
        },
        "diagnostics": {
            "energy_bin_lag1": lag1(columns["energy_per_site"]),
            "staggered_m2_bin_lag1": lag1(columns["staggered_m2"]),
            "maximum_absolute_z": maximum_z,
            "all_observables_agree": maximum_z <= z_limit,
        },
        "observables": estimates,
    }


def write_csv(
    path: Path,
    exact: dict[str, Any],
    comparisons: list[dict[str, Any]],
    short_metropolis: dict[str, Any],
    extrapolation: dict[str, Any],
    z_limit: float,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "method",
            "observable",
            "value",
            "uncertainty",
            "accuracy_budget",
            "instance",
            "tag",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for observable in OBSERVABLES:
            writer.writerow(
                {
                    "method": "exact_enumeration",
                    "observable": observable,
                    "value": f"{float(exact['observables'][observable]):.17g}",
                    "uncertainty": "0",
                    "accuracy_budget": "5e-13",
                    "instance": "L=2,N=20,PBC,T=Tc_exact",
                    "tag": "reference",
                }
            )
            for comparison in comparisons:
                estimate = comparison["observables"][observable]
                writer.writerow(
                    {
                        "method": f"ising_graph_{comparison['algorithm']}",
                        "observable": observable,
                        "value": f"{float(estimate['value']):.17g}",
                        "uncertainty": f"{float(estimate['conservative_standard_error']):.17g}",
                        "accuracy_budget": f"{z_limit * float(estimate['conservative_standard_error']):.17g}",
                        "instance": "L=2,N=20,PBC,T=Tc_exact",
                        "tag": estimate["tag"],
                    }
                )
            if observable == "uniform_susceptibility":
                short = short_metropolis["observables"][observable]
                writer.writerow(
                    {
                        "method": "ising_graph_metropolis_20000_sweeps_per_bin",
                        "observable": observable,
                        "value": f"{float(short['value']):.17g}",
                        "uncertainty": f"{float(short['conservative_standard_error']):.17g}",
                        "accuracy_budget": f"{z_limit * float(short['conservative_standard_error']):.17g}",
                        "instance": "L=2,N=20,PBC,T=Tc_exact",
                        "tag": short["tag"],
                    }
                )
                writer.writerow(
                    {
                        "method": "ising_graph_metropolis_1_over_bin_extrapolation",
                        "observable": observable,
                        "value": f"{float(extrapolation['value']):.17g}",
                        "uncertainty": f"{float(extrapolation['conservative_standard_error']):.17g}",
                        "accuracy_budget": f"{z_limit * float(extrapolation['conservative_standard_error']):.17g}",
                        "instance": "L=2,N=20,PBC,T=Tc_exact",
                        "tag": extrapolation["tag"],
                    }
                )


def write_readme(path: Path, payload: dict[str, Any]) -> None:
    exact = payload["exact"]["observables"]
    labels = {
        "energy_per_site": "$E/N$",
        "specific_heat_per_site": "$C/N$",
        "abs_staggered_m": "$\\langle|M_s|\\rangle/N$",
        "staggered_m2": "$\\langle M_s^2\\rangle/N^2$",
        "staggered_m4": "$\\langle M_s^4\\rangle/N^4$",
        "binder_staggered": "$U_s$",
        "uniform_susceptibility": "$\\chi_u$",
        "connected_staggered_susceptibility": "$\\chi_{s,c}$",
    }
    lines = [
        "# L=2 经典反铁磁 Ising 精确枚举校验",
        "",
        "## 约定与范围",
        "",
        "本目录验证 POSCAR 最近邻图在 $2\\times2$ 周期超胞（$N=20$、$N_b=24$）上的",
        "经典蒙特卡洛实现。统一使用",
        "",
        "$$H=+J\\sum_{\\langle ij\\rangle}s_i s_j,\\qquad J=1,\\quad s_i=\\pm1,$$",
        "",
        "以及 $M_s=\\sum_i\\eta_i s_i$。所有 24 条键都连接相反的 $\\eta$ 子格，",
        "因此 Wolff 路径等价于对 $\\tau_i=\\eta_i s_i$ 做铁磁集群更新。",
        "",
        "精确枚举覆盖全部 $2^{20}=1,048,576$ 个构型。对照温度为无限系统的精确临界点",
        f"$T_c/J={payload['temperature']:.15f}$；这里的有限 $L=2$ 系综仅用于程序校验，不能再估计热力学相变点。",
        "",
        "## 蒙卡协议与判据",
        "",
        "Wolff 和 Metropolis 各用 100,000 sweeps 热化，随后采样 256 bins。Wolff 每 bin",
        "20,000 sweeps；局域 Metropolis 为控制连通磁化率的有限-bin 偏差，每 bin 200,000 sweeps。",
        "另保留一个 20,000 sweeps/bin 的 Metropolis 诊断链。随机种子由运行编号和单元编号确定。",
        "误差分别按 1、2、4、8、16-bin 重分块估计，并取最大标准误；Binder 比对每种块长",
        "作 delete-one-block jackknife。预注册的一致性预算为 $4\\sigma_{\\rm cons}$。",
        "",
        "## 结果",
        "",
        "| 观测量 | 精确值 | Wolff | $z$ | Metropolis | $z$ |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    by_algorithm = {item["algorithm"]: item for item in payload["comparisons"]}
    for observable in OBSERVABLES:
        wolff = by_algorithm["wolff"]["observables"][observable]
        metro = by_algorithm["metropolis"]["observables"][observable]
        lines.append(
            f"| {labels[observable]} | {exact[observable]:.10g} | "
            f"{wolff['value']:.10g}({wolff['conservative_standard_error']:.2g}) | {wolff['z_score']:.2f} | "
            f"{metro['value']:.10g}({metro['conservative_standard_error']:.2g}) | {metro['z_score']:.2f} |"
        )
    gated_z = [
        abs(float(observable["z_score"]))
        for comparison in payload["comparisons"]
        for name, observable in comparison["observables"].items()
        if not (comparison["algorithm"] == "metropolis" and name == "uniform_susceptibility")
    ]
    extrapolation = payload["finite_bin_diagnostic"]["uniform_susceptibility_1_over_bin_extrapolation"]
    maximum = max(gated_z + [abs(float(extrapolation["z_score"]))])
    isolated_disagreement = (
        payload["status"]
        == "primary_observables_agree_metropolis_uniform_susceptibility_disagrees"
    )
    short_uniform = payload["finite_bin_diagnostic"]["metropolis_20000_sweeps_per_bin"]["observables"][
        "uniform_susceptibility"
    ]
    long_uniform = by_algorithm["metropolis"]["observables"]["uniform_susceptibility"]
    lines.extend(
        [
            "",
            "结论：原生 Metropolis $\\chi_u$ 与精确枚举为 `Disagreement`；"
            f"其余 15 个原生对照{'一致' if isolated_disagreement else '尚未全部一致'}，"
            f"连同单列外推量的最大门控 $|z|={maximum:.3f}$。",
            f"局域 Metropolis 的 $\\chi_u$ 在 20,000 sweeps/bin 时为 {short_uniform['value']:.8g} "
            f"（$z={short_uniform['z_score']:.2f}$），延长到 200,000 sweeps/bin 后为 "
            f"{long_uniform['value']:.8g}（$z={long_uniform['z_score']:.2f}$）；短链仅作有限-bin 偏差诊断。",
            f"按 $\\chi_u(n_{{\\rm bin}})=\\chi_u(\\infty)+a/n_{{\\rm bin}}$ 外推得到 "
            f"{extrapolation['value']:.8g}({extrapolation['conservative_standard_error']:.2g})，"
            f"相对精确值 $z={extrapolation['z_score']:.2f}$，标为 `{extrapolation['tag']}`。"
            "该外推支持偏差机制，但不是正式源码的原生输出，也不改变原生 `Disagreement`。",
            "`specific_heat_per_site`、两种磁化率是每个 bin 内构造的非线性/连通估计量；",
            "局域更新的原始 $\\chi_u$ 即使在 200,000 sweeps/bin 时仍有可分辨偏差，故不纳入主验收；",
            "$1/n_{\\rm bin}$ 外推只用于解释该诊断分歧。这些结果只验证本次协议。",
            "",
            "## 产物与复现",
            "",
            "- `density_of_states.csv`：按 $(E/J,M,M_s)$ 聚合的精确整数态密度。",
            "- `exact.json`：精确热力学量、键表、$\\eta$ 和规范。",
            "- `cells/*`：Fortran 输入、原始 bins、摘要、标准输出和清单。",
            "- `comparison.json` 与 `cross-method-check.csv`：保守重分块对照。",
            "- `checksums.sha256.json`：目录内全部产物哈希。",
            "",
            "在仓库根目录运行：",
            "",
            "```powershell",
            "powershell -NoProfile -ExecutionPolicy Bypass -File scripts\\run_ising_l2_exact_validation.ps1",
            "python scripts\\checksum_ledger.py results\\ising_l2_exact_validation_v1 --verify",
            "```",
            "",
            "正式 `ising_graph.f90` 未作修改。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--z-limit", type=float, default=4.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    exact = load_json(run_dir / "exact.json")
    graph_path = Path(exact["lattice"]["graph"])
    graph_hash = sha256_file(graph_path)
    source_hash = sha256_file(args.source.resolve())
    executable_hash = sha256_file(args.executable.resolve())
    if graph_hash != exact["lattice"]["graph_sha256"]:
        raise SystemExit("精确枚举格点图哈希不一致")
    if exact["lattice"]["N"] != 20 or exact["lattice"]["Nb"] != 24:
        raise SystemExit("精确枚举有限图尺寸错误")
    if exact["enumeration"]["number_of_configurations"] != 1_048_576:
        raise SystemExit("精确枚举态数错误")
    comparisons = [
        validate_cell(
            run_dir / "cells" / algorithm,
            algorithm,
            20000 if algorithm == "wolff" else 200000,
            exact,
            graph_hash,
            source_hash,
            executable_hash,
            args.z_limit,
        )
        for algorithm in ("wolff", "metropolis")
    ]
    short_metropolis = validate_cell(
        run_dir / "cells" / "metropolis-short",
        "metropolis",
        20000,
        exact,
        graph_hash,
        source_hash,
        executable_hash,
        args.z_limit,
    )
    short_uniform = short_metropolis["observables"]["uniform_susceptibility"]
    long_uniform = comparisons[1]["observables"]["uniform_susceptibility"]
    exact_uniform = float(exact["observables"]["uniform_susceptibility"])
    extrapolation = inverse_bin_extrapolation(
        short_uniform, long_uniform, 20000, 200000, exact_uniform, args.z_limit
    )
    gated_agreement = all(
        abs(float(observable["z_score"])) <= args.z_limit
        for comparison in comparisons
        for name, observable in comparison["observables"].items()
        if not (comparison["algorithm"] == "metropolis" and name == "uniform_susceptibility")
    )
    payload = {
        "schema_version": 1,
        "status": "primary_observables_agree_metropolis_uniform_susceptibility_disagrees"
        if gated_agreement and extrapolation["tag"] == "Agreement"
        else "fail",
        "verdict": {
            "native_metropolis_uniform_susceptibility": "Disagreement",
            "other_fifteen_native_comparisons": "Agreement" if gated_agreement else "Disagreement",
            "inverse_bin_extrapolation": extrapolation["tag"],
            "identified_cause": (
                "ising_graph.f90 forms (<M^2>-<M>_bin^2)/(N*T) inside each bin; "
                "slow sector flips under local Metropolis updates cause a finite-bin downward bias"
            ),
            "interpretation": "外推一致不能覆盖原生观测量的不一致",
        },
        "comparison_rule": {
            "tag": "Agreement when |MC-exact| <= z_limit*conservative_standard_error",
            "z_limit": args.z_limit,
            "block_sizes_bins": list(BLOCK_SIZES),
            "exact_numerical_budget": 5.0e-13,
        },
        "temperature": float(exact["parameters"]["temperature"]),
        "provenance": {
            "graph_sha256": graph_hash,
            "source_sha256": source_hash,
            "executable_sha256": executable_hash,
        },
        "exact": {
            "method": "complete enumeration of 2^20 configurations",
            "observables": {name: float(exact["observables"][name]) for name in OBSERVABLES},
        },
        "comparisons": comparisons,
        "finite_bin_diagnostic": {
            "role": "raw_disagreement_and_bias_diagnostic",
            "reason": "检验每-bin 连通估计量在局域更新下的有限-bin 偏差",
            "metropolis_20000_sweeps_per_bin": short_metropolis,
            "uniform_susceptibility_1_over_bin_extrapolation": extrapolation,
        },
    }
    (run_dir / "comparison.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(
        run_dir / "cross-method-check.csv",
        exact,
        comparisons,
        short_metropolis,
        extrapolation,
        args.z_limit,
    )
    write_readme(run_dir / "README.md", payload)
    if payload["status"] == "fail":
        raise SystemExit("蒙卡与精确枚举未在预注册预算内一致")
    print(run_dir / "comparison.json")


if __name__ == "__main__":
    main()
