#!/usr/bin/env python3
"""按 POSCAR 原坐标绘制周期最近邻图。"""

from __future__ import annotations

import argparse
from pathlib import Path

from lattice_from_poscar import parse_poscar, row_vector_times_matrix


def read_edges(path: Path) -> tuple[int, list[tuple[int, int, int, int]]]:
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    basis_count, edge_count, _ = (int(value) for value in lines[0])
    edges = []
    for fields in lines[1 : edge_count + 1]:
        site1, site2, dx, dy = (int(fields[index]) for index in range(4))
        edges.append((site1, site2, dx, dy))
    return basis_count, edges


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("poscar", type=Path)
    parser.add_argument("graph", type=Path)
    parser.add_argument("--cells", type=int, default=3)
    parser.add_argument("--output", type=Path, default=Path("results/poscar_lattice.svg"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    structure = parse_poscar(args.poscar)
    basis_count, edges = read_edges(args.graph)
    if basis_count != len(structure.atoms):
        raise SystemExit("POSCAR 格点数与图文件不一致")
    degree = {basis: 0 for basis in range(1, basis_count + 1)}
    for site1, site2, _, _ in edges:
        degree[site1] += 1
        degree[site2] += 1

    a1, a2 = structure.lattice[:2]
    local = {
        index: row_vector_times_matrix(atom.fractional, structure.lattice)[:2]
        for index, atom in enumerate(structure.atoms, start=1)
    }

    def position(cell_x: int, cell_y: int, basis: int) -> tuple[float, float]:
        point = local[basis]
        return point[0] + cell_x * a1[0] + cell_y * a2[0], point[1] + cell_x * a1[1] + cell_y * a2[1]

    cells = range(-args.cells, args.cells + 1)
    points = [position(x, y, basis) for x in cells for y in cells for basis in local]
    x_min, x_max = min(point[0] for point in points), max(point[0] for point in points)
    y_min, y_max = min(point[1] for point in points), max(point[1] for point in points)
    width, height = 1000, 760
    margin = 60.0
    scale = min((width - 2 * margin) / (x_max - x_min), (height - 2 * margin) / (y_max - y_min))

    def map_point(point: tuple[float, float]) -> tuple[float, float]:
        return margin + scale * (point[0] - x_min), height - margin - scale * (point[1] - y_min)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="500" y="30" text-anchor="middle" font-family="Arial" font-size="22">Direct rendering of POSCAR (2).vasp</text>',
    ]

    for cell_x in cells:
        for cell_y in cells:
            for site1, site2, dx, dy in edges:
                point1 = map_point(position(cell_x, cell_y, site1))
                point2 = map_point(position(cell_x + dx, cell_y + dy, site2))
                parts.append(
                    f'<line x1="{point1[0]:.2f}" y1="{point1[1]:.2f}" x2="{point2[0]:.2f}" y2="{point2[1]:.2f}" '
                    'stroke="#2458c6" stroke-width="2.4"/>'
                )

    for cell_x in cells:
        for cell_y in cells:
            for basis in local:
                x, y = map_point(position(cell_x, cell_y, basis))
                fill = "#aeb9cc" if degree[basis] == 3 else "#f5df48"
                radius = 5.8 if degree[basis] == 3 else 4.8
                parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.1f}" fill="{fill}" stroke="#173bb0" stroke-width="1.8"/>')

    offsets = {1: (13, 20), 2: (12, -10), 3: (12, -8), 4: (-25, -9), 5: (12, -9)}
    for basis in local:
        x, y = map_point(position(0, 0, basis))
        dx, dy = offsets.get(basis, (12, -9))
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="11" fill="none" stroke="#e0302d" stroke-width="2.4"/>')
        parts.append(f'<text x="{x + dx:.2f}" y="{y + dy:.2f}" font-family="Arial" font-size="18" fill="#c51f1a">{basis}</text>')

    origin = map_point(position(0, 0, 1))
    cell_corners = [
        origin,
        map_point((local[1][0] + a1[0], local[1][1] + a1[1])),
        map_point((local[1][0] + a1[0] + a2[0], local[1][1] + a1[1] + a2[1])),
        map_point((local[1][0] + a2[0], local[1][1] + a2[1])),
    ]
    polygon = " ".join(f"{x:.2f},{y:.2f}" for x, y in cell_corners)
    parts.append(f'<polygon points="{polygon}" fill="none" stroke="#e0302d" stroke-width="2" stroke-dasharray="8 5"/>')
    parts.append('<circle cx="80" cy="720" r="6" fill="#aeb9cc" stroke="#173bb0" stroke-width="1.5"/><text x="94" y="726" font-family="Arial" font-size="16">degree 3: sites 1, 2</text>')
    parts.append('<circle cx="300" cy="720" r="5" fill="#f5df48" stroke="#173bb0" stroke-width="1.5"/><text x="314" y="726" font-family="Arial" font-size="16">degree 2: sites 3, 4, 5</text>')
    parts.append("</svg>")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
