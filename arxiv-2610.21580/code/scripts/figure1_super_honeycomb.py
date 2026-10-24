#!/usr/bin/env python3
"""生成超蜂窝晶格与两类平带的论文 Figure 1。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7.0
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["xtick.major.width"] = 0.7
plt.rcParams["ytick.major.width"] = 0.7
plt.rcParams["xtick.major.size"] = 3.0
plt.rcParams["ytick.major.size"] = 3.0
plt.rcParams["legend.frameon"] = False


ROOT = Path(__file__).resolve().parents[1]
POSCAR = ROOT / "POSCAR (2).vasp"
GRAPH = ROOT / "graphs" / "super_honeycomb_poscar.dat"
TB_CSV = ROOT / "results" / "poscar2" / "bands_nn" / "band_path.csv"
TB_JSON = ROOT / "results" / "poscar2" / "bands_nn" / "flat_band_report.json"
LSWT_CSV = ROOT / "results" / "lswt_super_honeycomb_v1" / "magnon_path.csv"
LSWT_JSON = ROOT / "results" / "lswt_super_honeycomb_v1" / "summary.json"
OUTPUT = ROOT / "figures" / "figure1_super_honeycomb"

COLORS = {
    "A": "#2F6690",
    "B": "#E09F3E",
    "bond": "#AAB4BE",
    "band": "#4F6D8A",
    "band_light": "#8396A8",
    "goldstone": "#087E8B",
    "flat": "#B23A48",
    "grid": "#D9DEE3",
    "text": "#20252A",
}


def read_poscar(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """读取直接坐标并转为笛卡尔坐标。"""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    scale = float(lines[1])
    lattice = np.asarray([[float(value) for value in line.split()] for line in lines[2:5]]) * scale
    counts = [int(value) for value in lines[6].split()]
    coordinate_line = 7
    if lines[coordinate_line].lower().startswith("s"):
        coordinate_line += 1
    mode = lines[coordinate_line].lower()
    count = sum(counts)
    coordinates = np.asarray(
        [[float(value) for value in line.split()[:3]] for line in lines[coordinate_line + 1 : coordinate_line + 1 + count]]
    )
    if mode.startswith("d"):
        coordinates = coordinates @ lattice
    else:
        coordinates *= scale
    return lattice, coordinates


def read_graph(path: Path) -> tuple[list[tuple[int, int, int, int]], np.ndarray]:
    """读取周期图的键和二分子格标记。"""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    basis_count, edge_count, dimension = (int(value) for value in lines[0].split())
    if basis_count != 5 or dimension != 3:
        raise ValueError("Figure 1 需要五格点二维周期图")
    edges = []
    for line in lines[1 : edge_count + 1]:
        site1, site2, dx, dy, _, _ = (int(value) for value in line.split())
        edges.append((site1 - 1, site2 - 1, dx, dy))
    eta = np.asarray([int(value) for value in lines[edge_count + 1].split()])
    if tuple(np.count_nonzero(eta == value) for value in (1, -1)) != (2, 3):
        raise ValueError("二分子格计数不是 2:3")
    return edges, eta


def read_csv(path: Path) -> dict[str, np.ndarray]:
    """读取数值 CSV。"""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"空数据文件：{path}")
    return {name: np.asarray([float(row[name]) for row in rows]) for name in rows[0]}


def special_ticks(data: dict[str, np.ndarray]) -> tuple[list[float], list[str]]:
    """从路径坐标定位高对称点。"""
    targets = [(0.0, 0.0), (np.pi, 0.0), (2.0 * np.pi / 3.0, 2.0 * np.pi / 3.0), (0.0, 0.0)]
    distances = data["path_distance"]
    points = np.column_stack((data["kx"], data["ky"]))
    indices: list[int] = []
    start = 0
    for target_index, target in enumerate(targets):
        residual = np.linalg.norm(points - np.asarray(target), axis=1)
        if target_index == len(targets) - 1:
            index = len(points) - 1
        else:
            candidates = np.arange(start, len(points))
            index = int(candidates[np.argmin(residual[candidates])])
        indices.append(index)
        start = min(index + 1, len(points) - 1)
    return [float(distances[index]) for index in indices], [r"$\Gamma$", "M", "K", r"$\Gamma$"]


def style_band_axis(ax: plt.Axes, data: dict[str, np.ndarray], ylabel: str) -> None:
    """统一能带坐标轴。"""
    ticks, labels = special_ticks(data)
    for position in ticks:
        ax.axvline(position, color=COLORS["grid"], linewidth=0.65, zorder=0)
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_xticks(ticks, labels)
    ax.set_ylabel(ylabel, labelpad=2.5)
    ax.tick_params(direction="out", pad=2)
    ax.spines["left"].set_color("#555C62")
    ax.spines["bottom"].set_color("#555C62")


def draw_lattice(ax: plt.Axes, lattice: np.ndarray, basis: np.ndarray, edges, eta: np.ndarray) -> None:
    """绘制 POSCAR 最近邻周期图。"""
    a1, a2 = lattice[0, :2], lattice[1, :2]
    local = basis[:, :2]
    centre = 0.5 * (a1 + a2)
    xlim = (centre[0] - 10.6, centre[0] + 10.6)
    ylim = (centre[1] - 3.35, centre[1] + 3.35)

    for ix in range(-4, 5):
        for iy in range(-3, 4):
            origin = ix * a1 + iy * a2
            for site1, site2, dx, dy in edges:
                point1 = origin + local[site1]
                point2 = origin + dx * a1 + dy * a2 + local[site2]
                if (
                    max(point1[0], point2[0]) >= xlim[0] - 0.5
                    and min(point1[0], point2[0]) <= xlim[1] + 0.5
                    and max(point1[1], point2[1]) >= ylim[0] - 0.5
                    and min(point1[1], point2[1]) <= ylim[1] + 0.5
                ):
                    ax.plot(
                        [point1[0], point2[0]],
                        [point1[1], point2[1]],
                        color=COLORS["bond"],
                        linewidth=0.85,
                        solid_capstyle="round",
                        zorder=1,
                    )

    for ix in range(-4, 5):
        for iy in range(-3, 4):
            points = ix * a1 + iy * a2 + local
            visible = (
                (points[:, 0] >= xlim[0] - 0.3)
                & (points[:, 0] <= xlim[1] + 0.3)
                & (points[:, 1] >= ylim[0] - 0.3)
                & (points[:, 1] <= ylim[1] + 0.3)
            )
            for value, color in ((1, COLORS["A"]), (-1, COLORS["B"])):
                mask = visible & (eta == value)
                ax.scatter(
                    points[mask, 0],
                    points[mask, 1],
                    s=32,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.65,
                    zorder=3,
                )

    cell = np.asarray([[0.0, 0.0], a1, a1 + a2, a2, [0.0, 0.0]])
    ax.plot(cell[:, 0], cell[:, 1], color=COLORS["text"], linewidth=1.0, linestyle=(0, (3, 2)), zorder=4)
    for site, point in enumerate(local, start=1):
        text_color = "white" if eta[site - 1] == 1 else COLORS["text"]
        ax.text(
            point[0],
            point[1],
            str(site),
            ha="center",
            va="center",
            fontsize=5.6,
            fontweight="bold",
            color=text_color,
            zorder=5,
        )

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=5.4, markerfacecolor=COLORS["A"], markeredgecolor="white"),
        Line2D([], [], marker="o", linestyle="none", markersize=5.4, markerfacecolor=COLORS["B"], markeredgecolor="white"),
    ]
    ax.legend(
        handles,
        [r"$\mathcal{A}$: vertices (2/cell, $z_A=3$)", r"$\mathcal{B}$: edge centres (3/cell, $z_B=2$)"],
        loc="lower right",
        bbox_to_anchor=(0.995, 0.02),
        handletextpad=0.45,
        labelspacing=0.45,
        fontsize=6.7,
    )
    ax.text(
        0.018,
        0.10,
        r"sublattice imbalance  $N_A:N_B=2:3$",
        transform=ax.transAxes,
        fontsize=7.2,
        color=COLORS["text"],
        bbox={"boxstyle": "round,pad=0.27", "facecolor": "white", "edgecolor": "#C7CDD2", "linewidth": 0.65},
        zorder=7,
    )
    ax.text(0.012, 0.965, "a", transform=ax.transAxes, fontsize=9.0, fontweight="bold", va="top")
    ax.text(0.050, 0.965, "Five-site edge-centred honeycomb lattice", transform=ax.transAxes, fontsize=8.0, va="top")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_tight_binding(ax: plt.Axes, data: dict[str, np.ndarray], report: dict) -> None:
    """绘制单粒子最近邻能带。"""
    x = data["path_distance"]
    for name in ("band_1", "band_2", "band_4", "band_5"):
        ax.plot(x, data[name], color=COLORS["band"], linewidth=1.05, zorder=2)
    ax.plot(x, data["band_3"], color=COLORS["flat"], linewidth=2.2, zorder=4)
    flat = float(report["flat_energies"][0]["energy"])
    ax.annotate(
        r"flat band  $E/t=0$",
        xy=(0.72 * x[-1], flat),
        xytext=(0.55 * x[-1], 0.62),
        color=COLORS["flat"],
        fontsize=6.7,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["flat"], "linewidth": 0.75},
        ha="center",
    )
    style_band_axis(ax, data, r"Energy  $E/t$")
    ax.set_ylim(-2.62, 2.62)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_title("Nearest-neighbour hopping", loc="left", fontsize=8.0, pad=4)
    ax.text(-0.14, 1.035, "b", transform=ax.transAxes, fontsize=9.0, fontweight="bold", va="bottom")


def draw_lswt(ax: plt.Axes, data: dict[str, np.ndarray], report: dict) -> None:
    """绘制 LSWT 磁振子能带。"""
    x = data["path_distance"]
    quiet = ["omega_sigma_min_minus", "omega_sigma_max_plus", "omega_sigma_min_plus"]
    for index, name in enumerate(quiet):
        color = COLORS["band_light"] if index != 1 else COLORS["band"]
        ax.plot(x, data[name], color=color, linewidth=1.05, zorder=2)
    ax.plot(x, data["omega_sigma_max_minus"], color=COLORS["goldstone"], linewidth=1.75, zorder=3)
    ax.plot(x, data["omega_flat"], color=COLORS["flat"], linewidth=2.2, zorder=5)

    flat = float(report["bands"]["flat_energy_over_J"])
    exponent = float(report["goldstone"]["fitted_dynamic_exponent"])
    ax.annotate(
        rf"flat magnon  $\omega/J={flat:g}$",
        xy=(0.70 * x[-1], flat),
        xytext=(0.69 * x[-1], 1.30),
        color=COLORS["flat"],
        fontsize=6.7,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["flat"], "linewidth": 0.75},
        ha="center",
    )
    goldstone_x = 0.085 * x[-1]
    goldstone_y = float(np.interp(goldstone_x, x, data["omega_sigma_max_minus"]))
    ax.annotate(
        "type-B Goldstone\n" + rf"$z={exponent:.3f}$",
        xy=(goldstone_x, goldstone_y),
        xytext=(0.24 * x[-1], 0.23),
        color=COLORS["goldstone"],
        fontsize=6.7,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["goldstone"], "linewidth": 0.75},
        ha="center",
        va="center",
    )
    style_band_axis(ax, data, r"Magnon energy  $\omega/J$")
    ax.set_ylim(0.0, 1.58)
    ax.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax.set_title(r"Heisenberg antiferromagnet (LSWT, $S=1/2$)", loc="left", fontsize=8.0, pad=4)
    ax.text(-0.14, 1.035, "c", transform=ax.transAxes, fontsize=9.0, fontweight="bold", va="bottom")


def validate_inputs(tb_data: dict[str, np.ndarray], tb_report: dict, lswt_report: dict) -> None:
    """核对图中核心数值。"""
    if len(tb_report.get("flat_energies", [])) != 1:
        raise ValueError("单粒子平带计数异常")
    if abs(float(tb_report["flat_energies"][0]["energy"])) > 1.0e-10:
        raise ValueError("单粒子平带不在零能")
    if np.ptp(tb_data["band_3"]) > 1.0e-10:
        raise ValueError("band_3 不是严格平带")
    if int(lswt_report["bands"]["flat_band_count"]) != 1:
        raise ValueError("LSWT 平带计数异常")
    if abs(float(lswt_report["bands"]["flat_energy_over_J"]) - 1.0) > 1.0e-10:
        raise ValueError("LSWT 平带能量异常")
    if abs(float(lswt_report["goldstone"]["fitted_dynamic_exponent"]) - 2.0) > 0.01:
        raise ValueError("Goldstone 指数偏离 z=2")


def main() -> int:
    lattice, basis = read_poscar(POSCAR)
    edges, eta = read_graph(GRAPH)
    tb_data = read_csv(TB_CSV)
    lswt_data = read_csv(LSWT_CSV)
    tb_report = json.loads(TB_JSON.read_text(encoding="utf-8"))
    lswt_report = json.loads(LSWT_JSON.read_text(encoding="utf-8"))
    validate_inputs(tb_data, tb_report, lswt_report)

    width_mm, height_mm = 183.0, 140.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.06, 1.0],
        left=0.075,
        right=0.985,
        bottom=0.095,
        top=0.985,
        hspace=0.28,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])

    draw_lattice(ax_a, lattice, basis, edges, eta)
    draw_tight_binding(ax_b, tb_data, tb_report)
    draw_lswt(ax_c, lswt_data, lswt_report)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".svg"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".tiff"), dpi=600, facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Figure 1 已写入：{OUTPUT.parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
