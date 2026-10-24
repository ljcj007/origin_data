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


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 9.3
plt.rcParams["mathtext.fontset"] = "dejavusans"
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["xtick.major.width"] = 0.7
plt.rcParams["ytick.major.width"] = 0.7
plt.rcParams["xtick.major.size"] = 3.0
plt.rcParams["ytick.major.size"] = 3.0
plt.rcParams["legend.frameon"] = False


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"
POSCAR = DATA / "POSCAR (2).vasp"
GRAPH = DATA / "super_honeycomb_poscar.dat"
TB_CSV = DATA / "band_path.csv"
TB_JSON = DATA / "flat_band_report.json"
LSWT_CSV = DATA / "magnon_path.csv"
LSWT_JSON = DATA / "summary.json"
OUTPUT = ROOT / "figures" / "figure1_super_honeycomb"

COLORS = {
    "A": "#336990",
    "B": "#D99A38",
    "bond": "#BAC3CC",
    "band": "#516F8B",
    "band_light": "#899CAB",
    "goldstone": "#137F88",
    "flat": "#AF3C50",
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
    ax.set_ylabel(ylabel, labelpad=3.0)
    ax.tick_params(direction="out", pad=2)
    ax.spines["left"].set_color("#555C62")
    ax.spines["bottom"].set_color("#555C62")


def draw_lattice(ax: plt.Axes, lattice: np.ndarray, basis: np.ndarray, edges, eta: np.ndarray) -> None:
    """保留真实连接，突出五格点原胞。"""
    scale = np.linalg.norm(lattice[0, :2])
    a1, a2 = lattice[:2, :2] / scale
    local = basis[:, :2] / scale
    xlim, ylim = (-1.75, 2.25), (-0.37, 1.23)

    def visible(point):
        return xlim[0] + 0.06 < point[0] < xlim[1] - 0.06 and ylim[0] + 0.06 < point[1] < ylim[1] - 0.06

    for ix in range(-3, 4):
        for iy in range(-2, 3):
            origin = ix * a1 + iy * a2
            for site1, site2, dx, dy in edges:
                p = origin + local[site1]
                q = origin + dx * a1 + dy * a2 + local[site2]
                if visible(p) and visible(q):
                    ax.plot(*np.array([p, q]).T, color=COLORS["bond"], lw=1.05, zorder=1)
            for index, point in enumerate(origin + local):
                if visible(point):
                    color = COLORS["A"] if eta[index] == 1 else COLORS["B"]
                    ax.scatter(*point, s=32, c=color, edgecolors="white", linewidths=0.5, zorder=3)

    cell = np.asarray([[0.0, 0.0], a1, a1 + a2, a2, [0.0, 0.0]])
    ax.plot(*cell.T, color=COLORS["text"], lw=0.9, ls=(0, (3.5, 2.5)), zorder=2)
    for index, point in enumerate(local):
        color = COLORS["A"] if eta[index] == 1 else COLORS["B"]
        ax.scatter(*point, s=69, c=color, edgecolors="white", linewidths=0.65, zorder=4)
        ax.text(*point, str(index + 1), ha="center", va="center", fontsize=8.2,
                color="white" if eta[index] == 1 else "#252A30", fontweight="bold", zorder=5)
    ax.set(xlim=xlim, ylim=ylim, aspect="equal")
    ax.axis("off")


def draw_lattice_key(ax: plt.Axes) -> None:
    """图例与晶格分开排布。"""
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    for y, color, title, subtitle in [
        (0.86, COLORS["A"], r"$\mathcal{A}$  Hub sites", r"2 per cell, $z_A=3$"),
        (0.51, COLORS["B"], r"$\mathcal{B}$  Edge sites", r"3 per cell, $z_B=2$"),
    ]:
        ax.scatter(0.035, y, s=43, color=color, edgecolors="white", linewidths=0.5)
        ax.text(0.13, y, title, va="center", fontsize=9.3, color=COLORS["text"])
        ax.text(0.13, y - 0.145, subtitle, va="center", fontsize=8.6, color="#555C62")
    ax.plot([0.005, 0.09], [0.145, 0.145], color=COLORS["text"], lw=0.9, ls=(0, (3.5, 2.5)))
    ax.text(0.13, 0.145, "Primitive cell", fontsize=8.8, va="center")


def draw_tight_binding(ax: plt.Axes, data: dict[str, np.ndarray], report: dict) -> None:
    """绘制单粒子最近邻能带。"""
    x = data["path_distance"]
    for name in ("band_1", "band_2", "band_4", "band_5"):
        ax.plot(x, data[name], color=COLORS["band"], linewidth=1.05, zorder=2)
    ax.plot(x, data["band_3"], color=COLORS["flat"], linewidth=2.2, zorder=4)
    flat = float(report["flat_energies"][0]["energy"])
    ax.annotate(
        r"flat band  $E/t=0$",
        xy=(0.69 * x[-1], flat),
        xytext=(0.52 * x[-1], 0.65),
        color=COLORS["flat"],
        fontsize=9.0,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["flat"], "linewidth": 0.75},
        ha="center",
    )
    style_band_axis(ax, data, r"Energy  $E/t$")
    ax.set_ylim(-2.62, 2.62)
    ax.set_yticks([-2, -1, 0, 1, 2])


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
        xy=(0.82 * x[-1], flat),
        xytext=(0.56 * x[-1], 1.40),
        color=COLORS["flat"],
        fontsize=9.0,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["flat"], "linewidth": 0.75},
        ha="center",
    )
    goldstone_x = 0.085 * x[-1]
    goldstone_y = float(np.interp(goldstone_x, x, data["omega_sigma_max_minus"]))
    ax.annotate(
        "type-B Goldstone\n" + rf"$z={exponent:.3f}$",
        xy=(goldstone_x, goldstone_y),
        xytext=(0.56 * x[-1], 0.19),
        color=COLORS["goldstone"],
        fontsize=8.8,
        arrowprops={"arrowstyle": "-|>", "color": COLORS["goldstone"], "linewidth": 0.75},
        ha="center",
        va="center",
    )
    style_band_axis(ax, data, r"$\omega/J$")
    ax.set_ylim(0.0, 1.58)
    ax.set_yticks([0.0, 0.5, 1.0, 1.5])


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

    width_mm, height_mm = 174.0, 106.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    ax_a = fig.add_axes([0.015, 0.563, 0.68, 0.353])
    ax_key = fig.add_axes([0.720, 0.586, 0.27, 0.302])
    ax_b = fig.add_axes([0.082, 0.086, 0.395, 0.398])
    ax_c = fig.add_axes([0.586, 0.086, 0.395, 0.398])

    draw_lattice(ax_a, lattice, basis, edges, eta)
    draw_lattice_key(ax_key)
    draw_tight_binding(ax_b, tb_data, tb_report)
    draw_lswt(ax_c, lswt_data, lswt_report)
    for x, y, letter, title in [
        (0.022, 0.955, "(a)", "Edge-centered honeycomb lattice"),
        (0.022, 0.520, "(b)", "Single-particle bands"),
        (0.527, 0.520, "(c)", r"LSWT magnons ($S=1/2$)"),
    ]:
        fig.text(x, y, letter, fontsize=10.0, fontweight="bold", va="center")
        fig.text(x + 0.06, y, title, fontsize=9.7, va="center")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".svg"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Figure 1 已写入：{OUTPUT.parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
