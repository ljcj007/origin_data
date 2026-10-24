#!/usr/bin/env python3
"""绘制 POSCAR 已确认的 super-honeycomb 晶格。"""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "poscar2" / "lattice_structure.svg"

A1 = (1.5, -math.sqrt(3.0) / 2.0)
A2 = (1.5, math.sqrt(3.0) / 2.0)
POSITIONS = {
    1: (0.0, 0.0),
    2: (1.0, 0.0),
    3: (0.5, 0.0),
    4: (-0.25, math.sqrt(3.0) / 4.0),
    5: (-0.25, -math.sqrt(3.0) / 4.0),
}
LABELS = {1: "A", 2: "B", 3: "C", 4: "D", 5: "E"}
EDGES = [
    (1, 3, 0, 0),
    (2, 3, 0, 0),
    (1, 4, 0, 0),
    (2, 4, 1, 0),
    (1, 5, 0, 0),
    (2, 5, 0, 1),
]


def position(cell_x: int, cell_y: int, basis: int) -> tuple[float, float]:
    local = POSITIONS[basis]
    return (
        cell_x * A1[0] + cell_y * A2[0] + local[0],
        cell_x * A1[1] + cell_y * A2[1] + local[1],
    )


def main() -> int:
    cells = range(-3, 4)
    raw_points = [position(x, y, basis) for x in cells for y in cells for basis in POSITIONS]
    x_values = [point[0] for point in raw_points]
    y_values = [point[1] for point in raw_points]
    width, height = 1000, 680
    margin = 55.0
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    scale = min((width - 2 * margin) / (x_max - x_min), (height - 2 * margin) / (y_max - y_min))

    def map_point(point: tuple[float, float]) -> tuple[float, float]:
        return (
            margin + scale * (point[0] - x_min),
            height - margin - scale * (point[1] - y_min),
        )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="500" y="30" text-anchor="middle" font-family="Arial" font-size="22">POSCAR-confirmed edge-centered honeycomb lattice</text>',
    ]

    for cell_x in cells:
        for cell_y in cells:
            for site1, site2, dx, dy in EDGES:
                point1 = map_point(position(cell_x, cell_y, site1))
                point2 = map_point(position(cell_x + dx, cell_y + dy, site2))
                parts.append(
                    f'<line x1="{point1[0]:.2f}" y1="{point1[1]:.2f}" '
                    f'x2="{point2[0]:.2f}" y2="{point2[1]:.2f}" stroke="#356db0" stroke-width="2.2"/>'
                )

    for cell_x in cells:
        for cell_y in cells:
            for basis in POSITIONS:
                x, y = map_point(position(cell_x, cell_y, basis))
                color = "#2f77b4" if basis <= 2 else "#f3c74f"
                parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5.3" fill="{color}" stroke="#202020" stroke-width="0.8"/>')

    central = [map_point(position(0, 0, basis)) for basis in (4, 1, 5, 3, 2)]
    polygon = " ".join(f"{x:.2f},{y:.2f}" for x, y in central)
    parts.append(f'<polygon points="{polygon}" fill="none" stroke="#e33b32" stroke-width="2.5" stroke-dasharray="8 5"/>')
    for basis in POSITIONS:
        x, y = map_point(position(0, 0, basis))
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="10" fill="none" stroke="#e33b32" stroke-width="2"/>')
        parts.append(f'<text x="{x + 12:.2f}" y="{y - 10:.2f}" font-family="Arial" font-size="18" fill="#b51f1a">{LABELS[basis]}</text>')

    origin = map_point(position(0, 0, 1))
    for vector, label in ((A1, "a₁"), (A2, "a₂")):
        endpoint = map_point((vector[0], vector[1]))
        parts.append(
            f'<line x1="{origin[0]:.2f}" y1="{origin[1]:.2f}" x2="{endpoint[0]:.2f}" y2="{endpoint[1]:.2f}" '
            'stroke="#222" stroke-width="2" marker-end="url(#arrow)"/>'
        )
        parts.append(f'<text x="{endpoint[0] + 8:.2f}" y="{endpoint[1]:.2f}" font-family="Arial" font-size="17">{label}</text>')

    parts.insert(
        2,
        '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#222"/></marker></defs>',
    )
    parts.append('<circle cx="70" cy="625" r="6" fill="#2f77b4"/><text x="84" y="631" font-family="Arial" font-size="16">honeycomb vertices (2/cell)</text>')
    parts.append('<circle cx="340" cy="625" r="6" fill="#f3c74f" stroke="#202020"/><text x="354" y="631" font-family="Arial" font-size="16">edge centers (3/cell)</text>')
    parts.append("</svg>")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
