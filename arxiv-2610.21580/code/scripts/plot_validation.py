"""汇总 beta 扫描并生成无第三方绘图库的 SVG。"""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BETAS = [2, 4, 8, 16]


def read_summary(path: Path) -> dict[str, tuple[float, float]]:
    data: dict[str, tuple[float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 3:
            try:
                data[parts[0]] = (float(parts[1]), float(parts[2]))
            except ValueError:
                pass
    return data


def scale(value: float, low: float, high: float, start: float, size: float) -> float:
    return start + (value - low) * size / (high - low)


def panel(
    x0: int,
    y0: int,
    width: int,
    height: int,
    values: list[float],
    errors: list[float],
    exact: float,
    title: str,
    ylabel: str,
) -> list[str]:
    ymin = min([exact] + [v - e for v, e in zip(values, errors)])
    ymax = max([exact] + [v + e for v, e in zip(values, errors)])
    margin = 0.08 * (ymax - ymin)
    ymin -= margin
    ymax += margin
    left, top = x0 + 70, y0 + 35
    plot_w, plot_h = width - 95, height - 90

    def px(i: int) -> float:
        return left + i * plot_w / (len(BETAS) - 1)

    def py(v: float) -> float:
        return top + plot_h - scale(v, ymin, ymax, 0, plot_h)

    lines = [
        f'<text x="{x0 + width/2}" y="{y0 + 20}" class="title">{title}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{py(exact):.2f}" x2="{left + plot_w}" y2="{py(exact):.2f}" class="exact"/>',
        f'<text x="{left + plot_w - 4}" y="{py(exact) - 6:.2f}" class="exact-label">ED</text>',
        f'<text x="{x0 + 16}" y="{top + plot_h/2}" class="ylabel" transform="rotate(-90 {x0 + 16} {top + plot_h/2})">{ylabel}</text>',
        f'<text x="{left + plot_w/2}" y="{y0 + height - 12}" class="label">beta</text>',
    ]

    for tick in range(5):
        value = ymin + tick * (ymax - ymin) / 4
        y = py(value)
        lines.append(f'<line x1="{left - 5}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" class="axis"/>')
        lines.append(f'<text x="{left - 9}" y="{y + 4:.2f}" class="tick" text-anchor="end">{value:.4f}</text>')

    points: list[str] = []
    for i, (beta, value, error) in enumerate(zip(BETAS, values, errors)):
        x, y = px(i), py(value)
        y1, y2 = py(value - error), py(value + error)
        lines.extend(
            [
                f'<line x1="{x:.2f}" y1="{y1:.2f}" x2="{x:.2f}" y2="{y2:.2f}" class="error"/>',
                f'<line x1="{x - 4:.2f}" y1="{y1:.2f}" x2="{x + 4:.2f}" y2="{y1:.2f}" class="error"/>',
                f'<line x1="{x - 4:.2f}" y1="{y2:.2f}" x2="{x + 4:.2f}" y2="{y2:.2f}" class="error"/>',
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" class="point"/>',
                f'<text x="{x:.2f}" y="{top + plot_h + 19}" class="tick">{beta}</text>',
            ]
        )
        points.append(f"{x:.2f},{y:.2f}")
    lines.append(f'<polyline points="{" ".join(points)}" class="guide"/>')
    return lines


def main() -> None:
    validation = ROOT / "results" / "validation"
    rows = []
    for beta in BETAS:
        summary = read_summary(validation / f"beta_{beta}" / "results_summary.txt")
        energy, energy_err = summary["energy_per_site"]
        m2, m2_err = summary["staggered_m2"]
        rows.append((beta, energy, energy_err, m2, m2_err))

    ed = json.loads((ROOT / "results" / "ed_L2.json").read_text(encoding="utf-8"))
    csv_path = validation / "validation_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["beta", "energy", "energy_stderr", "staggered_m2", "staggered_m2_stderr"])
        writer.writerows(rows)

    width, height = 920, 430
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        """<style>
        text{font-family:Arial,sans-serif;fill:#20242b}.title{font-size:17px;font-weight:700;text-anchor:middle}
        .label{font-size:14px;text-anchor:middle}.ylabel{font-size:14px;text-anchor:middle}.tick{font-size:12px;text-anchor:middle}
        .axis{stroke:#20242b;stroke-width:1.2}.error{stroke:#1769aa;stroke-width:1.4}.point{fill:#1769aa}
        .guide{fill:none;stroke:#1769aa;stroke-width:1.3;opacity:.55}.exact{stroke:#d64b32;stroke-width:1.6;stroke-dasharray:6 4}
        .exact-label{fill:#d64b32;font-size:12px;text-anchor:end}.note{font-size:13px}
        </style>""",
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    svg += panel(10, 15, 450, 365, [r[1] for r in rows], [r[2] for r in rows], ed["ground_state_energy_per_site"], "Energy convergence", "E/N")
    svg += panel(460, 15, 450, 365, [r[3] for r in rows], [r[4] for r in rows], ed["ground_state_staggered_m2"], "Staggered magnetization", "m_s^2")
    svg.append('<text x="460" y="414" class="note" text-anchor="middle">Lx=Ly=2, N=16, J=1; error bars are bin standard errors</text>')
    svg.append("</svg>")
    (validation / "beta_convergence.svg").write_text("\n".join(svg), encoding="utf-8")
    print(csv_path)
    print(validation / "beta_convergence.svg")


if __name__ == "__main__":
    main()
