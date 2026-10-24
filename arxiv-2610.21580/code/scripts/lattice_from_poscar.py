#!/usr/bin/env python3
"""分析 POSCAR 的距离壳层，并导出周期格点图。"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


Vector = tuple[float, float, float]
Shift = tuple[int, int, int]


@dataclass(frozen=True)
class Atom:
    original_index: int
    species: str
    fractional: Vector


@dataclass(frozen=True)
class Structure:
    title: str
    lattice: tuple[Vector, Vector, Vector]
    atoms: tuple[Atom, ...]


@dataclass(frozen=True)
class Edge:
    site1: int
    site2: int
    shift: Shift
    distance: float


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def determinant(matrix: Sequence[Sequence[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def inverse(matrix: Sequence[Sequence[float]]) -> tuple[Vector, Vector, Vector]:
    a, b, c = matrix
    det = determinant(matrix)
    if abs(det) < 1.0e-14:
        raise ValueError("晶格矩阵不可逆")
    return (
        (
            (b[1] * c[2] - b[2] * c[1]) / det,
            (a[2] * c[1] - a[1] * c[2]) / det,
            (a[1] * b[2] - a[2] * b[1]) / det,
        ),
        (
            (b[2] * c[0] - b[0] * c[2]) / det,
            (a[0] * c[2] - a[2] * c[0]) / det,
            (a[2] * b[0] - a[0] * b[2]) / det,
        ),
        (
            (b[0] * c[1] - b[1] * c[0]) / det,
            (a[1] * c[0] - a[0] * c[1]) / det,
            (a[0] * b[1] - a[1] * b[0]) / det,
        ),
    )


def row_vector_times_matrix(vector: Sequence[float], matrix: Sequence[Sequence[float]]) -> Vector:
    return tuple(sum(vector[i] * matrix[i][j] for i in range(3)) for j in range(3))  # type: ignore[return-value]


def parse_poscar(path: Path) -> Structure:
    lines = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if len(lines) < 8:
        raise ValueError("POSCAR 内容不足")

    title = lines[0]
    scale = float(lines[1].split()[0])
    raw_lattice = tuple(tuple(float(x) for x in lines[i].split()[:3]) for i in range(2, 5))
    if any(len(v) != 3 for v in raw_lattice):
        raise ValueError("晶格矢量必须有三个分量")

    if scale < 0.0:
        volume = abs(determinant(raw_lattice))
        if volume <= 0.0:
            raise ValueError("负缩放因子需要非零晶胞体积")
        factor = (-scale / volume) ** (1.0 / 3.0)
    else:
        factor = scale
    lattice = tuple(tuple(factor * x for x in row) for row in raw_lattice)

    cursor = 5
    tokens = lines[cursor].split()
    has_symbols = not all(token.lstrip("+-").isdigit() for token in tokens)
    if has_symbols:
        symbols = tokens
        cursor += 1
        counts = [int(x) for x in lines[cursor].split()]
    else:
        counts = [int(x) for x in tokens]
        symbols = [f"X{i + 1}" for i in range(len(counts))]
    if len(symbols) != len(counts):
        raise ValueError("元素名数量与原子计数不一致")
    cursor += 1

    if lines[cursor].lower().startswith("s"):
        cursor += 1
    coordinate_mode = lines[cursor].lower()
    cursor += 1
    if not coordinate_mode.startswith(("d", "c", "k")):
        raise ValueError("坐标类型必须是 Direct 或 Cartesian")

    labels = [symbol for symbol, count in zip(symbols, counts) for _ in range(count)]
    if len(lines) < cursor + len(labels):
        raise ValueError("原子坐标数量不足")

    inv_lattice = inverse(lattice)
    atoms: list[Atom] = []
    for index, species in enumerate(labels, start=1):
        values = tuple(float(x) for x in lines[cursor + index - 1].split()[:3])
        if len(values) != 3:
            raise ValueError(f"第 {index} 个原子坐标格式错误")
        if coordinate_mode.startswith("d"):
            fractional = values
        else:
            cartesian = tuple(factor * x for x in values)
            fractional = row_vector_times_matrix(cartesian, inv_lattice)
        fractional = tuple(x - math.floor(x) for x in fractional)
        atoms.append(Atom(index, species, fractional))

    return Structure(title, lattice, tuple(atoms))


def metric_lattice(structure: Structure, metric: str) -> tuple[Vector, Vector, Vector]:
    if metric == "cartesian":
        return structure.lattice
    if metric == "normalized":
        normalized: list[Vector] = []
        for vector in structure.lattice:
            length = math.sqrt(dot(vector, vector))
            if length <= 0.0:
                raise ValueError("存在零长度晶格矢量")
            normalized.append(tuple(x / length for x in vector))
        return tuple(normalized)  # type: ignore[return-value]
    if metric == "fractional":
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    raise ValueError(f"未知度量：{metric}")


def parse_periodic(text: str) -> tuple[bool, bool, bool]:
    clean = text.lower().replace(",", "").replace(" ", "")
    if clean == "none":
        return False, False, False
    if any(axis not in "xyz" for axis in clean) or len(set(clean)) != len(clean):
        raise argparse.ArgumentTypeError("周期方向应写成 xy、xyz 或 none")
    return tuple(axis in clean for axis in "xyz")  # type: ignore[return-value]


def canonical_edge(site1: int, site2: int, shift: Shift) -> tuple[int, int, Shift]:
    forward = (site1, site2, shift)
    reverse = (site2, site1, tuple(-x for x in shift))
    return min(forward, reverse)


def enumerate_edges(
    atoms: Sequence[Atom],
    lattice: Sequence[Sequence[float]],
    periodic: Sequence[bool],
    image_range: int,
) -> list[Edge]:
    image_axes: list[Iterable[int]] = [
        range(-image_range, image_range + 1) if enabled else (0,) for enabled in periodic
    ]
    unique: dict[tuple[int, int, Shift], Edge] = {}

    for site1, atom1 in enumerate(atoms, start=1):
        for site2, atom2 in enumerate(atoms, start=1):
            for shift in itertools.product(*image_axes):
                shift_tuple = tuple(int(x) for x in shift)
                if site1 == site2 and shift_tuple == (0, 0, 0):
                    continue
                key = canonical_edge(site1, site2, shift_tuple)
                if key in unique:
                    continue
                delta = tuple(atom2.fractional[i] + shift_tuple[i] - atom1.fractional[i] for i in range(3))
                cartesian = row_vector_times_matrix(delta, lattice)
                distance = math.sqrt(dot(cartesian, cartesian))
                if distance > 1.0e-12:
                    unique[key] = Edge(key[0], key[1], key[2], distance)

    return sorted(unique.values(), key=lambda edge: (edge.distance, edge.site1, edge.site2, edge.shift))


def group_shells(edges: Sequence[Edge], absolute_tolerance: float, relative_tolerance: float) -> list[list[Edge]]:
    shells: list[list[Edge]] = []
    for edge in edges:
        if not shells:
            shells.append([edge])
            continue
        reference = sum(item.distance for item in shells[-1]) / len(shells[-1])
        tolerance = max(absolute_tolerance, relative_tolerance * max(reference, edge.distance))
        if abs(edge.distance - reference) <= tolerance:
            shells[-1].append(edge)
        else:
            shells.append([edge])
    return shells


def degree_statistics(site_count: int, edges: Sequence[Edge]) -> tuple[int, int, float, list[int]]:
    degrees = [0] * site_count
    for edge in edges:
        degrees[edge.site1 - 1] += 1
        degrees[edge.site2 - 1] += 1
    return min(degrees), max(degrees), sum(degrees) / site_count, degrees


def solve_bipartition(
    site_count: int, edges: Sequence[Edge], periodic: Sequence[bool]
) -> tuple[bool, list[int], Shift | None]:
    active_axes = [index for index, enabled in enumerate(periodic) if enabled]
    for phase_values in itertools.product((0, 1), repeat=len(active_axes)):
        phase = [0, 0, 0]
        for axis, value in zip(active_axes, phase_values):
            phase[axis] = value
        colors: list[int | None] = [None] * site_count
        consistent = True
        changed = True
        while changed and consistent:
            changed = False
            for edge in edges:
                i = edge.site1 - 1
                j = edge.site2 - 1
                translation_parity = sum(phase[k] * edge.shift[k] for k in range(3)) % 2
                required = 1 ^ translation_parity
                if colors[i] is None and colors[j] is None:
                    colors[i] = 0
                    colors[j] = required
                    changed = True
                elif colors[i] is None:
                    colors[i] = colors[j] ^ required  # type: ignore[operator]
                    changed = True
                elif colors[j] is None:
                    colors[j] = colors[i] ^ required
                    changed = True
                elif colors[i] ^ colors[j] != required:
                    consistent = False
                    break
        if consistent:
            finished = [0 if color is None else color for color in colors]
            eta = [1 if color == 0 else -1 for color in finished]
            return True, eta, tuple(phase)  # type: ignore[return-value]
    return False, [0] * site_count, None


def select_atoms(structure: Structure, species: Sequence[str] | None) -> tuple[Atom, ...]:
    if not species:
        return structure.atoms
    requested = set(species)
    available = {atom.species for atom in structure.atoms}
    missing = sorted(requested - available)
    if missing:
        raise ValueError(f"POSCAR 中没有这些元素：{', '.join(missing)}")
    return tuple(atom for atom in structure.atoms if atom.species in requested)


def edge_record(edge: Edge) -> dict[str, object]:
    return {
        "site1": edge.site1,
        "site2": edge.site2,
        "shift": list(edge.shift),
        "distance": edge.distance,
    }


def analyze_metric(
    structure: Structure,
    atoms: Sequence[Atom],
    metric: str,
    periodic: Sequence[bool],
    image_range: int,
    maximum_shells: int,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> dict[str, object]:
    edges = enumerate_edges(atoms, metric_lattice(structure, metric), periodic, image_range)
    shells = group_shells(edges, absolute_tolerance, relative_tolerance)[:maximum_shells]
    records: list[dict[str, object]] = []
    cumulative: list[Edge] = []
    for number, shell in enumerate(shells, start=1):
        cumulative.extend(shell)
        degree_min, degree_max, degree_mean, degrees = degree_statistics(len(atoms), cumulative)
        bipartite, eta, phase = solve_bipartition(len(atoms), cumulative, periodic)
        records.append(
            {
                "shell": number,
                "distance_mean": sum(edge.distance for edge in shell) / len(shell),
                "distance_min": min(edge.distance for edge in shell),
                "distance_max": max(edge.distance for edge in shell),
                "edge_count_in_shell": len(shell),
                "edge_count_cumulative": len(cumulative),
                "degree_cumulative": {
                    "minimum": degree_min,
                    "maximum": degree_max,
                    "mean": degree_mean,
                    "by_site": degrees,
                },
                "bipartite_cumulative": bipartite,
                "eta": eta,
                "translation_phase": list(phase) if phase is not None else None,
                "edges": [edge_record(edge) for edge in shell],
            }
        )
    return {"metric": metric, "shells": records}


def print_report(report: dict[str, object]) -> None:
    print(f"结构：{report['title']}")
    print(f"选中格点：{report['site_count']} / {report['atom_count']}")
    print(f"元素：{', '.join(report['selected_species'])}")
    print(f"周期方向：{report['periodic']}")
    for metric_report in report["metrics"]:  # type: ignore[index]
        print(f"\n度量：{metric_report['metric']}")
        print("壳层       距离        本壳键数   累计配位范围   二分图")
        for shell in metric_report["shells"]:
            degree = shell["degree_cumulative"]
            print(
                f"{shell['shell']:>4}  {shell['distance_mean']:>11.7f}"
                f"  {shell['edge_count_in_shell']:>10}"
                f"   {degree['minimum']}-{degree['maximum']}"
                f"            {'是' if shell['bipartite_cumulative'] else '否'}"
            )


def create_report(args: argparse.Namespace) -> tuple[Structure, tuple[Atom, ...], dict[str, object]]:
    structure = parse_poscar(args.poscar)
    atoms = select_atoms(structure, args.species)
    if not atoms:
        raise ValueError("没有选中任何格点")
    metrics = [
        analyze_metric(
            structure,
            atoms,
            metric,
            args.periodic,
            args.image_range,
            args.max_shells,
            args.absolute_tolerance,
            args.relative_tolerance,
        )
        for metric in args.metric
    ]
    periodic_text = "".join(axis for axis, enabled in zip("xyz", args.periodic) if enabled) or "none"
    report: dict[str, object] = {
        "format": "poscar-lattice-analysis-v1",
        "source": str(args.poscar.resolve()),
        "title": structure.title,
        "atom_count": len(structure.atoms),
        "site_count": len(atoms),
        "selected_species": sorted({atom.species for atom in atoms}),
        "periodic": periodic_text,
        "image_range": args.image_range,
        "atoms": [
            {
                "site": site,
                "original_index": atom.original_index,
                "species": atom.species,
                "fractional": list(atom.fractional),
            }
            for site, atom in enumerate(atoms, start=1)
        ],
        "metrics": metrics,
    }
    return structure, atoms, report


def run_analyze(args: argparse.Namespace) -> None:
    _, _, report = create_report(args)
    print_report(report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\n完整报告：{args.output}")


def run_build(args: argparse.Namespace) -> None:
    if len(args.metric) != 1:
        raise ValueError("导出连边时只能选择一种度量")
    _, atoms, report = create_report(args)
    metric_report = report["metrics"][0]  # type: ignore[index]
    available_shells = metric_report["shells"]
    if args.shells < 1 or args.shells > len(available_shells):
        raise ValueError("请求的壳层数量超出分析范围")

    chosen: list[dict[str, object]] = available_shells[: args.shells]
    edges: list[dict[str, object]] = []
    for shell in chosen:
        for edge in shell["edges"]:
            edges.append({**edge, "bond_type": shell["shell"]})
    bipartite, eta, phase = solve_bipartition(
        len(atoms),
        [
            Edge(item["site1"], item["site2"], tuple(item["shift"]), item["distance"])  # type: ignore[arg-type]
            for item in edges
        ],
        args.periodic,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(atoms)} {len(edges)} 3\n")
        for edge in edges:
            dx, dy, dz = edge["shift"]
            handle.write(
                f"{edge['site1']} {edge['site2']} {dx} {dy} {dz} {edge['bond_type']}\n"
            )
        handle.write(" ".join(str(value) for value in eta) + "\n")
        handle.write(" ".join(str(value) for value in (phase or (0, 0, 0))) + "\n")

    companion = args.output.with_suffix(args.output.suffix + ".json")
    build_record = {
        **report,
        "format": "periodic-lattice-graph-v1",
        "chosen_metric": args.metric[0],
        "chosen_shells": args.shells,
        "edge_count": len(edges),
        "bipartite": bipartite,
        "eta": eta,
        "translation_phase": list(phase) if phase is not None else None,
        "graph_file": str(args.output.resolve()),
    }
    companion.write_text(json.dumps(build_record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print_report(report)
    print(f"\n连边文件：{args.output}")
    print(f"构建记录：{companion}")
    print(f"二分图：{'是' if bipartite else '否'}")


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("poscar", type=Path, help="POSCAR 文件")
    parser.add_argument("--periodic", type=parse_periodic, default=parse_periodic("xy"), help="周期方向")
    parser.add_argument("--species", nargs="*", help="作为格点的元素；默认全部")
    parser.add_argument(
        "--metric",
        nargs="+",
        choices=("cartesian", "normalized", "fractional"),
        default=("cartesian", "normalized"),
        help="距离度量",
    )
    parser.add_argument("--image-range", type=int, default=2, help="搜索的周期像范围")
    parser.add_argument("--max-shells", type=int, default=8, help="报告的壳层数")
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-7, help="壳层绝对容差")
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-5, help="壳层相对容差")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser("analyze", help="只报告距离壳层")
    add_common_arguments(analyze)
    analyze.add_argument("--output", type=Path, default=Path("results/lattice_analysis.json"))
    analyze.set_defaults(function=run_analyze)

    build = subparsers.add_parser("build", help="按已确认壳层导出连边")
    add_common_arguments(build)
    build.set_defaults(metric=("cartesian",))
    build.add_argument("--shells", type=int, required=True, help="累计采用的壳层数")
    build.add_argument("--output", type=Path, default=Path("lattice_graph.dat"))
    build.set_defaults(function=run_build)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.image_range < 1 or args.max_shells < 1:
        parser.error("image-range 和 max-shells 必须为正整数")
    try:
        args.function(args)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
