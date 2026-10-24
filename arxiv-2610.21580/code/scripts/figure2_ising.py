#!/usr/bin/env python3
"""生成 Wolff Ising 临界温度验证图。"""

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
plt.rcParams["font.size"] = 7.0
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["xtick.major.width"] = 0.7
plt.rcParams["ytick.major.width"] = 0.7
plt.rcParams["xtick.major.size"] = 3.0
plt.rcParams["ytick.major.size"] = 3.0


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "results" / "ising_tc_wolff_v1"
AUDIT_DIR = ROOT / "results" / "ising_tc_wolff_audit_v1"
DATA_CSV = AUDIT_DIR / "audit_observables.csv"
SUMMARY_JSON = AUDIT_DIR / "audit_summary.json"
OUTPUT = ROOT / "figures" / "figure2_ising"
SIZES = (16, 24, 32, 48, 64)

COLORS = {
    16: "#A9C6DA",
    24: "#77A9C6",
    32: "#4B88AF",
    48: "#28678F",
    64: "#123F5A",
}
MARKERS = {16: "o", 24: "s", 32: "D", 48: "^", 64: "v"}
EXACT_COLOR = "#B23A48"


def read_rows(path: Path) -> list[dict[str, float | int | str]]:
    """读取正式扫描汇总表。"""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    parsed = []
    for row in rows:
        parsed.append(
            {
                "cell_id": row["cell_id"],
                "temperature": float(row["temperature"]),
                "L": int(row["L"]),
                "binder": float(row["binder"]),
                "binder_error": float(row["binder_error_conservative"]),
            }
        )
    return parsed


def validate(rows: list[dict[str, float | int | str]], summary: dict) -> None:
    """核对数据完整性和误差定义所需条件。"""
    if len(rows) != 35 or int(summary["cell_count"]) != 35:
        raise ValueError("Wolff 正式扫描不是 35 个单元")
    if summary.get("status") != "pass_with_provenance_boundary":
        raise ValueError("Wolff 严格审计未通过")
    sizes = sorted({int(row["L"]) for row in rows})
    temperatures = sorted({float(row["temperature"]) for row in rows})
    if sizes != [8, 12, 16, 24, 32, 48, 64] or len(temperatures) != 5:
        raise ValueError("尺寸或温度轴不完整")
    if any(float(row["binder_error"]) <= 0.0 for row in rows):
        raise ValueError("Binder 误差必须为正")
    exact_formula = 2.0 / np.arccosh(2.0 + np.sqrt(3.0))
    if abs(float(summary["exact_temperature"]) - exact_formula) > 1.0e-12:
        raise ValueError("精确临界温度不一致")


def selected_crossings(summary: dict) -> list[dict[str, float | int]]:
    """提取四组大尺寸相邻交点。"""
    selected = []
    allowed = {(16, 24), (24, 32), (32, 48), (48, 64)}
    for item in summary["conservative_crossings"]:
        pair = (int(item["L_small"]), int(item["L_large"]))
        if pair not in allowed:
            continue
        roots = item["roots"]
        if len(roots) != 1:
            raise ValueError(f"crossing 根数量异常：{pair}")
        root = roots[0]
        selected.append(
            {
                "L_small": pair[0],
                "L_large": pair[1],
                "temperature": float(root["temperature"]),
                "error": float(root["statistical_error"]),
                "inverse_size": 1.0 / np.sqrt(pair[0] * pair[1]),
            }
        )
    if len(selected) != 4:
        raise ValueError("缺少大尺寸 crossing")
    return sorted(selected, key=lambda item: float(item["inverse_size"]))


def style_axis(ax: plt.Axes) -> None:
    """统一坐标轴细节。"""
    ax.tick_params(direction="out", pad=2)
    ax.spines["left"].set_color("#555C62")
    ax.spines["bottom"].set_color("#555C62")


def draw_binder(ax: plt.Axes, rows: list[dict[str, float | int | str]], exact: float) -> None:
    """绘制 Binder 曲线。"""
    for size in SIZES:
        series = sorted((row for row in rows if int(row["L"]) == size), key=lambda row: float(row["temperature"]))
        x = np.asarray([float(row["temperature"]) for row in series])
        y = np.asarray([float(row["binder"]) for row in series])
        error = np.asarray([float(row["binder_error"]) for row in series])
        ax.errorbar(
            x,
            y,
            yerr=error,
            color=COLORS[size],
            marker=MARKERS[size],
            markersize=3.5,
            markerfacecolor="white" if size < 48 else COLORS[size],
            markeredgewidth=0.8,
            linewidth=1.15,
            elinewidth=0.7,
            capsize=1.8,
            capthick=0.7,
            zorder=3,
        )
        ax.text(1.01235, y[-1], rf"$L={size}$", color=COLORS[size], fontsize=6.2, va="center")

    ax.axvline(exact, color=EXACT_COLOR, linewidth=1.0, linestyle=(0, (3, 2)), zorder=1)
    ax.text(
        exact + 0.00022,
        0.491,
        rf"exact $T_c/J={exact:.6f}$",
        color=EXACT_COLOR,
        fontsize=6.2,
        rotation=90,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 0.7},
    )
    ax.set_xlim(0.9954, 1.0138)
    ax.set_ylim(0.475, 0.661)
    ax.set_xticks([0.996, 1.000, 1.004, 1.008, 1.012])
    ax.set_yticks([0.48, 0.52, 0.56, 0.60, 0.64])
    ax.set_ylabel(r"Staggered Binder cumulant  $U_4$")
    ax.set_xlabel(r"Temperature  $T/J$", labelpad=2)
    ax.text(-0.16, 1.015, "a", transform=ax.transAxes, fontsize=9.0, fontweight="bold", va="bottom")
    ax.set_title("Wolff finite-size crossings", loc="left", fontsize=8.0, pad=4)
    style_axis(ax)


def draw_crossings(ax: plt.Axes, crossings: list[dict[str, float | int]], exact: float) -> None:
    """绘制交点的有限尺寸趋近。"""
    x = np.asarray([float(item["inverse_size"]) for item in crossings])
    y = np.asarray([float(item["temperature"]) for item in crossings])
    error = np.asarray([float(item["error"]) for item in crossings])
    ax.plot(x, y, color="#8AA7B8", linewidth=0.8, linestyle=(0, (2, 2)), zorder=1)
    ax.errorbar(
        x,
        y,
        yerr=error,
        color="#28678F",
        marker="o",
        markersize=3.8,
        markerfacecolor="white",
        markeredgewidth=0.85,
        linewidth=0.0,
        elinewidth=0.8,
        capsize=2.2,
        capthick=0.8,
        zorder=3,
    )
    ax.axhline(exact, color=EXACT_COLOR, linewidth=1.0, linestyle=(0, (3, 2)), zorder=0)
    ax.text(0.0545, exact + 0.00006, "exact", color=EXACT_COLOR, fontsize=6.0, ha="right", va="bottom")

    offsets = [(0.0005, 0.00034), (0.0006, 0.00039), (0.0007, -0.00043), (-0.0005, -0.00048)]
    alignments = [("left", "bottom"), ("left", "bottom"), ("left", "top"), ("right", "top")]
    for item, offset, alignment in zip(crossings, offsets, alignments):
        label = f"{item['L_small']}–{item['L_large']}"
        ax.text(
            float(item["inverse_size"]) + offset[0],
            float(item["temperature"]) + offset[1],
            label,
            fontsize=5.7,
            color="#3D4A52",
            ha=alignment[0],
            va=alignment[1],
        )

    ax.annotate(
        r"$L\rightarrow\infty$",
        xy=(0.0135, 1.00215),
        xytext=(0.023, 1.00215),
        fontsize=6.2,
        color="#4B565D",
        arrowprops={"arrowstyle": "-|>", "linewidth": 0.7, "color": "#4B565D"},
        ha="center",
        va="center",
    )
    ax.set_xlim(0.012, 0.056)
    ax.set_ylim(1.0018, 1.00515)
    ax.set_xticks([0.02, 0.03, 0.04, 0.05])
    ax.set_yticks([1.002, 1.003, 1.004, 1.005])
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_xlabel(r"Inverse pair size  $1/\sqrt{L_1L_2}$", labelpad=2)
    ax.set_ylabel(r"Crossing  $T_\times/J$")
    ax.text(-0.16, 1.015, "b", transform=ax.transAxes, fontsize=9.0, fontweight="bold", va="bottom")
    ax.set_title("Large-size convergence", loc="left", fontsize=8.0, pad=4)
    style_axis(ax)


def main() -> int:
    rows = read_rows(DATA_CSV)
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
    validate(rows, summary)
    crossings = selected_crossings(summary)
    exact = float(summary["exact_temperature"])

    width_mm, height_mm = 89.0, 116.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=[2.05, 1.0],
        left=0.19,
        right=0.955,
        bottom=0.105,
        top=0.965,
        hspace=0.43,
    )
    draw_binder(fig.add_subplot(grid[0]), rows, exact)
    draw_crossings(fig.add_subplot(grid[1]), crossings, exact)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".svg"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".tiff"), dpi=600, facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Figure 2 已写入：{OUTPUT.parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
