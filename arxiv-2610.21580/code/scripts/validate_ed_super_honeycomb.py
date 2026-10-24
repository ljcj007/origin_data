#!/usr/bin/env python3
"""严格验收 2x2 super-honeycomb 海森堡 ED 结果。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import deque
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH = PROJECT_ROOT / "graphs" / "super_honeycomb_poscar.dat"
DEFAULT_RESULT = PROJECT_ROOT / "results" / "ed_super_honeycomb_v1" / "result.json"
EXPECTED_GRAPH_SHA256 = "58BE9C1BC747B0B17F45B31216E096E82CEBFB97874B3AAE5ADCC2DFAE1D65FF"
EXPECTED_LX = 2
EXPECTED_LY = 2
EXPECTED_NSITE = 20
EXPECTED_NBOND = 24
EXPECTED_NPLUS = 8
EXPECTED_NMINUS = 12
EXPECTED_SPIN = 2.0
EXPECTED_SPIN_SQUARED = 6.0
EXPECTED_FERRIMAGNETIC_M2 = 0.015
EXPECTED_TOL = 1.0e-11
EXPECTED_MAXITER = 1200
EXPECTED_SSE_SEED = 350837008
ENERGY_TOLERANCE = 1.0e-9
RESIDUAL_LIMIT = 1.0e-9
STATE_NORM_TOLERANCE = 1.0e-10
SPIN_SQUARED_TOLERANCE = 1.0e-8
MOMENT_TOLERANCE = 1.0e-10


class ValidationError(RuntimeError):
    """表示 ED 验收条件未满足。"""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def integer(value: Any, name: str) -> int:
    require(isinstance(value, int) and not isinstance(value, bool), f"{name} 必须为整数")
    return int(value)


def finite(value: Any, name: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool), f"{name} 必须为数值")
    result = float(value)
    require(math.isfinite(result), f"{name} 必须为有限数")
    return result


def close(actual: float, expected: float, tolerance: float, name: str) -> None:
    require(abs(actual - expected) <= tolerance, f"{name} 不一致：{actual}，期望 {expected}")


def graph_records(path: Path) -> tuple[int, list[tuple[int, int, int, int]], list[int], list[int]]:
    require(path.is_file(), f"缺少格点图：{path}")
    lines = [
        content.split()
        for raw in path.read_text(encoding="utf-8").splitlines()
        if (content := raw.split("#", 1)[0].strip())
    ]
    require(bool(lines) and len(lines[0]) == 3, "格点图首行格式错误")
    n_basis, n_edges, dimension = (int(value) for value in lines[0])
    require((n_basis, n_edges, dimension) == (5, 6, 3), "目标格点图头信息不符")
    require(len(lines) == n_edges + 3, "格点图记录数量错误")

    edges: list[tuple[int, int, int, int]] = []
    for tokens in lines[1 : n_edges + 1]:
        require(len(tokens) == 6, "格点图键记录格式错误")
        site1, site2, dx, dy, dz, bond_type = (int(value) for value in tokens)
        require(1 <= site1 <= n_basis and 1 <= site2 <= n_basis, "格点图含越界编号")
        require(dz == 0 and bond_type == 1, "目标 ED 仅接受 xy 面内的单一最近邻键")
        edges.append((site1 - 1, site2 - 1, dx, dy))

    eta = [int(value) for value in lines[n_edges + 1]]
    phase = [int(value) for value in lines[n_edges + 2]]
    require(len(eta) == n_basis and set(eta) <= {-1, 1}, "二分子格标记错误")
    require(len(phase) == 3 and set(phase) <= {0, 1}, "子格平移相位错误")
    return n_basis, edges, eta, phase


def inspect_graph(path: Path, lx: int = EXPECTED_LX, ly: int = EXPECTED_LY) -> dict[str, Any]:
    digest = sha256_file(path)
    require(digest == EXPECTED_GRAPH_SHA256, "目标格点图 SHA-256 不一致")
    n_basis, edge_defs, eta_basis, phase = graph_records(path)
    nsite = n_basis * lx * ly

    def site(x: int, y: int, local: int) -> int:
        return n_basis * ((x % lx) + (y % ly) * lx) + local

    eta = [0] * nsite
    bonds: list[tuple[int, int]] = []
    for y in range(ly):
        for x in range(lx):
            parity = (phase[0] * x + phase[1] * y) % 2
            sign = -1 if parity else 1
            for local in range(n_basis):
                eta[site(x, y, local)] = eta_basis[local] * sign
            for local1, local2, dx, dy in edge_defs:
                bonds.append((site(x, y, local1), site(x + dx, y + dy, local2)))

    require(all(i != j for i, j in bonds), "周期图含自键")
    undirected = {(min(i, j), max(i, j)) for i, j in bonds}
    require(len(undirected) == len(bonds), "周期图含重复无向键")
    require(all(eta[i] * eta[j] == -1 for i, j in bonds), "周期图不是二分图")

    adjacency = [[] for _ in range(nsite)]
    for i, j in bonds:
        adjacency[i].append(j)
        adjacency[j].append(i)
    reached = {0}
    queue = deque([0])
    while queue:
        for neighbor in adjacency[queue.popleft()]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    require(len(reached) == nsite, "周期图不连通")

    nplus = eta.count(1)
    nminus = eta.count(-1)
    require((nsite, len(bonds), nplus, nminus) == (20, 24, 8, 12), "目标簇计数错误")
    require(sorted(set(map(len, adjacency))) == [2, 3], "目标簇配位数错误")
    return {
        "path": str(path.resolve()),
        "sha256": digest,
        "lx": lx,
        "ly": ly,
        "nsite": nsite,
        "nbond": len(bonds),
        "n_plus": nplus,
        "n_minus": nminus,
        "degrees": sorted(set(map(len, adjacency))),
        "connected": True,
        "bipartite": True,
    }


def validate_sectors(result: dict[str, Any], nsite: int) -> dict[str, Any]:
    raw = result.get("sector_energies")
    require(isinstance(raw, list) and len(raw) == 4, "必须包含 Sz=0,1,2,3 四个扇区")
    sectors: dict[int, dict[str, Any]] = {}
    for entry in raw:
        require(isinstance(entry, dict), "扇区记录格式错误")
        magnetization = integer(entry.get("magnetization"), "magnetization")
        require(magnetization in range(4) and magnetization not in sectors, "Sz 扇区重复或越界")
        nup = integer(entry.get("nup"), f"Sz={magnetization} nup")
        random_seed = integer(entry.get("random_seed"), f"Sz={magnetization} random_seed")
        dimension = integer(entry.get("dimension"), f"Sz={magnetization} dimension")
        expected_dimension = math.comb(nsite, nsite // 2 + magnetization)
        require(nup == nsite // 2 + magnetization, f"Sz={magnetization} 的 nup 错误")
        require(
            random_seed == 20260903 + nup,
            f"Sz={magnetization} 的随机种子错误",
        )
        require(dimension == expected_dimension, f"Sz={magnetization} 的块维数错误")
        require(
            integer(entry.get("expected_dimension"), f"Sz={magnetization} expected_dimension")
            == expected_dimension,
            f"Sz={magnetization} 的预期块维数错误",
        )
        energy = finite(entry.get("energy"), f"Sz={magnetization} energy")
        residual = finite(entry.get("residual"), f"Sz={magnetization} residual")
        state_norm = finite(entry.get("state_norm"), f"Sz={magnetization} state_norm")
        require(0.0 <= residual <= RESIDUAL_LIMIT, f"Sz={magnetization} Lanczos 残差超限")
        close(state_norm, 1.0, STATE_NORM_TOLERANCE, f"Sz={magnetization} 态范数")
        sectors[magnetization] = {
            "nup": nup,
            "random_seed": random_seed,
            "dimension": dimension,
            "energy": energy,
            "residual": residual,
            "state_norm": state_norm,
        }
    require(set(sectors) == {0, 1, 2, 3}, "Sz 扇区集合不完整")

    energy0 = sectors[0]["energy"]
    spread = max(abs(sectors[sz]["energy"] - energy0) for sz in (0, 1, 2))
    gap = sectors[3]["energy"] - energy0
    require(spread <= ENERGY_TOLERANCE, "Sz=0,1,2 未在 1e-9 内简并")
    require(gap > 0.0, "Sz=3 相对基态的能隙必须为正")
    close(
        finite(result.get("ground_multiplet_energy_spread"), "ground_multiplet_energy_spread"),
        spread,
        1.0e-12,
        "多重态能量展宽",
    )
    close(
        finite(result.get("gap_to_next_spin_sector"), "gap_to_next_spin_sector"),
        gap,
        1.0e-12,
        "Sz=3 能隙",
    )
    close(finite(result.get("ground_state_energy"), "ground_state_energy"), energy0, 1.0e-12, "基态能量")
    close(
        finite(result.get("ground_state_energy_per_site"), "ground_state_energy_per_site"),
        energy0 / nsite,
        1.0e-12,
        "每站点基态能量",
    )
    residual_max = max(sector["residual"] for sector in sectors.values())
    close(
        finite(result.get("maximum_lanczos_residual"), "maximum_lanczos_residual"),
        residual_max,
        1.0e-15,
        "最大 Lanczos 残差",
    )
    return {
        "sectors": [{"magnetization": sz, **sectors[sz]} for sz in range(4)],
        "ground_multiplet_energy_spread": spread,
        "gap_to_sz3": gap,
        "maximum_lanczos_residual": residual_max,
    }


def validate_observables(result: dict[str, Any], nsite: int) -> dict[str, float]:
    ground_spin = finite(result.get("lieb_mattis_ground_spin"), "lieb_mattis_ground_spin")
    close(ground_spin, EXPECTED_SPIN, 0.0, "Lieb-Mattis 基态总自旋")
    declared_s2 = finite(result.get("lieb_mattis_spin_squared"), "lieb_mattis_spin_squared")
    close(declared_s2, EXPECTED_SPIN_SQUARED, 0.0, "Lieb-Mattis S(S+1)")
    measured_s2 = finite(result.get("uniform_spin_squared"), "uniform_spin_squared")
    close(measured_s2, EXPECTED_SPIN_SQUARED, SPIN_SQUARED_TOLERANCE, "<S_tot^2>")
    ferrimagnetic_m2 = finite(result.get("ferrimagnetic_m2"), "ferrimagnetic_m2")
    close(ferrimagnetic_m2, measured_s2 / nsite**2, 1.0e-12, "mf2 与 <S_tot^2>/N^2")
    close(ferrimagnetic_m2, EXPECTED_FERRIMAGNETIC_M2, MOMENT_TOLERANCE, "mf2 精确目标")
    staggered_m2 = finite(result.get("staggered_m2"), "staggered_m2")
    require(0.0 <= staggered_m2 <= 0.25 + 0.5 / nsite, "staggered_m2 超出自旋上界")
    return {
        "lieb_mattis_ground_spin": ground_spin,
        "uniform_spin_squared": measured_s2,
        "ferrimagnetic_m2": ferrimagnetic_m2,
        "staggered_m2": staggered_m2,
    }


def compare_sse_pilot(manifest: dict[str, Any], ed_result: dict[str, Any], graph_hash: str) -> dict[str, Any]:
    require(manifest.get("schema_version") == 1, "SSE 清单版本错误")
    require(manifest.get("status", {}).get("success") is True, "SSE pilot 未成功完成")
    params = manifest.get("params", {})
    settings = manifest.get("settings", {})
    require(integer(params.get("L"), "SSE L") == 2, "SSE pilot 不是同一 L=2 簇")
    close(finite(params.get("beta"), "SSE beta"), 8.0, 0.0, "SSE pilot beta")
    close(finite(params.get("beta_over_L2"), "SSE beta_over_L2"), 2.0, 0.0, "SSE beta/L^2")
    require(settings.get("model") == "spin_half_antiferromagnetic_heisenberg", "SSE 模型不一致")
    require(settings.get("boundary") == "periodic_xy", "SSE 边界条件不一致")
    require(settings.get("couplings") == [1.0], "SSE 耦合不一致")
    require(
        (settings.get("basis_sites"), settings.get("basis_plus"), settings.get("basis_minus"))
        == (5, 2, 3),
        "SSE 晶胞计数不一致",
    )
    provenance = manifest.get("provenance", {})
    execution = manifest.get("execution", {})
    require(str(provenance.get("graph_sha256", "")).upper() == graph_hash, "SSE 格点图哈希不一致")
    require(integer(execution.get("returncode"), "SSE returncode") == 0, "SSE 返回码非零")
    require(integer(execution.get("seed"), "SSE seed") == EXPECTED_SSE_SEED, "SSE 随机种子不一致")

    observables = manifest.get("observables", {})
    energy = finite(observables.get("energy_per_site", {}).get("value"), "SSE energy_per_site")
    energy_error = finite(observables.get("energy_per_site", {}).get("error"), "SSE energy error")
    require(energy_error > 0.0, "SSE 能量误差必须为正")
    ed_energy = finite(ed_result.get("ground_state_energy_per_site"), "ED energy_per_site")
    lower_bound_z = (energy - ed_energy) / energy_error
    require(lower_bound_z >= -5.0, "SSE 有限温能量显著低于 ED 基态能量")

    comparisons: dict[str, Any] = {
        "beta": 8.0,
        "energy_per_site": {
            "ed_zero_temperature": ed_energy,
            "sse_finite_temperature": energy,
            "sse_standard_error": energy_error,
            "difference_in_sse_sigma": lower_bound_z,
            "ground_energy_bound_consistent": True,
        },
    }
    for name in ("ferrimagnetic_m2", "staggered_m2"):
        value = finite(observables.get(name, {}).get("value"), f"SSE {name}")
        error = finite(observables.get(name, {}).get("error"), f"SSE {name} error")
        require(error > 0.0, f"SSE {name} 误差必须为正")
        target = finite(ed_result.get(name), f"ED {name}")
        comparisons[name] = {
            "ed_zero_temperature": target,
            "sse_finite_temperature": value,
            "sse_standard_error": error,
            "difference_in_sse_sigma": (value - target) / error,
        }
    comparisons["interpretation"] = "SSE 为 beta=8 的有限温结果；除基态能量下界外，差值仅报告而不作为零温一致性硬门。"
    return comparisons


def validate_result(
    result: dict[str, Any],
    graph_path: Path = DEFAULT_GRAPH,
    sse_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    require(result.get("schema_version") == 2, "ED 结果 schema_version 必须为 2")
    require(result.get("status") == "completed", "ED 结果状态不是 completed")
    require(result.get("hamiltonian") == "H = sum_<ij> S_i dot S_j, J=1", "ED 哈密顿量约定不一致")
    solver = result.get("solver", {})
    require(solver.get("method") == "scipy.sparse.linalg.eigsh", "ED 求解器记录不一致")
    require(solver.get("random_seed_rule") == "20260903+nup", "ED 随机种子规则不一致")
    close(finite(solver.get("tol"), "solver tol"), EXPECTED_TOL, 0.0, "Lanczos tol")
    require(
        integer(solver.get("maxiter"), "solver maxiter") == EXPECTED_MAXITER,
        "Lanczos maxiter 不一致",
    )

    graph = inspect_graph(graph_path)
    require(str(result.get("graph", "")).replace("\\", "/").endswith("graphs/super_honeycomb_poscar.dat"), "ED 图路径声明不一致")
    require(str(result.get("graph_sha256", "")).upper() == graph["sha256"], "ED 图哈希不一致")
    generator = PROJECT_ROOT / "scripts" / "ed_graph_heisenberg.py"
    require(str(result.get("source_sha256", "")).upper() == sha256_file(generator), "ED 源码哈希不一致")
    for key, expected in (
        ("lx", EXPECTED_LX),
        ("ly", EXPECTED_LY),
        ("nsite", EXPECTED_NSITE),
        ("nbond", EXPECTED_NBOND),
        ("n_plus", EXPECTED_NPLUS),
        ("n_minus", EXPECTED_NMINUS),
    ):
        require(integer(result.get(key), key) == expected, f"{key} 不一致")
        graph_key = {"nsite": "nsite", "nbond": "nbond", "n_plus": "n_plus", "n_minus": "n_minus"}.get(key)
        if graph_key:
            require(integer(result.get(key), key) == graph[graph_key], f"{key} 与格点图不一致")

    sectors = validate_sectors(result, graph["nsite"])
    observables = validate_observables(result, graph["nsite"])
    report: dict[str, Any] = {
        "schema_version": 1,
        "status": "passed",
        "graph": graph,
        "sectors": sectors,
        "observables": observables,
        "acceptance": {
            "multiplet_energy_tolerance": ENERGY_TOLERANCE,
            "lanczos_residual_limit": RESIDUAL_LIMIT,
            "state_norm_tolerance": STATE_NORM_TOLERANCE,
            "spin_squared_tolerance": SPIN_SQUARED_TOLERANCE,
            "ferrimagnetic_m2_tolerance": MOMENT_TOLERANCE,
        },
    }
    if sse_manifest is not None:
        report["sse_pilot_comparison"] = compare_sse_pilot(sse_manifest, result, graph["sha256"])
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--sse-manifest", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = json.loads(args.result.read_text(encoding="utf-8"))
        sse_manifest = (
            json.loads(args.sse_manifest.read_text(encoding="utf-8"))
            if args.sse_manifest is not None
            else None
        )
        report = validate_result(result, args.graph.resolve(), sse_manifest)
    except (OSError, ValueError, KeyError, ValidationError) as exc:
        raise SystemExit(f"ED 验证失败：{exc}") from exc

    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
