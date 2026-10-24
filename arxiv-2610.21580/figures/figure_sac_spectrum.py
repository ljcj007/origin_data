#!/usr/bin/env python3
"""绘制 QMC+SAC 动力学谱及 LSWT 参照线。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

local_packages = Path(__file__).resolve().parents[1] / ".python-packages"
if local_packages.is_dir():
    sys.path.insert(0, str(local_packages))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def cell_edges(values: np.ndarray) -> np.ndarray:
    if len(values) == 1:
        return np.asarray([values[0] - 0.5, values[0] + 0.5])
    middle = 0.5 * (values[1:] + values[:-1])
    return np.concatenate(([values[0] - (middle[0] - values[0])], middle, [values[-1] + (values[-1] - middle[-1])]))


def spectrum_grid(
    rows: list[dict[str, str]], channel: str, manifest: dict
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    chosen = [row for row in rows if row["channel"] == channel]
    q_indices = [int(point["q_index"]) for point in manifest["q_points"]]
    omega = np.asarray(sorted({float(row["omega"]) for row in chosen}))
    distance_values = []
    cumulative = 0.0
    previous: tuple[float, float] | None = None
    for point in manifest["q_points"]:
        current = (float(point["kx"]), float(point["ky"]))
        if previous is not None:
            cumulative += math.dist(previous, current)
        distance_values.append(cumulative)
        previous = current
    distance = np.asarray(distance_values)
    grid = np.zeros((len(omega), len(q_indices)))
    q_lookup = {value: index for index, value in enumerate(q_indices)}
    w_lookup = {round(value, 12): index for index, value in enumerate(omega)}
    present = {int(row["q_index"]) for row in chosen}
    for q_index in q_indices:
        if q_index not in present:
            grid[:, q_lookup[q_index]] = np.nan
    for row in chosen:
        grid[w_lookup[round(float(row["omega"]), 12)], q_lookup[int(row["q_index"])]] = float(row["S_q_omega"])
    return distance, omega, grid


def tick_positions(manifest: dict, rows: list[dict[str, str]]) -> tuple[list[float], list[str]]:
    distance = {}
    for row in rows:
        distance[int(row["q_index"])] = float(row["path_distance"])
    positions: list[float] = []
    labels: list[str] = []
    for point in manifest["q_points"]:
        if point["label"] != "path":
            positions.append(distance[point["q_index"]])
            labels.append(point["label"])
    return positions, labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectra", type=Path)
    parser.add_argument("channels_manifest", type=Path)
    parser.add_argument("lswt_path", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--channels", default="total,edge,cls")
    parser.add_argument("--omega-max", type=float, default=2.0)
    parser.add_argument("--title-suffix", default="pilot, L=6")
    args = parser.parse_args()

    rows = read_csv(args.spectra)
    manifest = json.loads(args.channels_manifest.read_text(encoding="utf-8"))
    lswt = read_csv(args.lswt_path)
    channels = [value.strip() for value in args.channels.split(",") if value.strip()]
    channel_titles = {
        "total": r"physical $S(\mathbf{q},\omega)$",
        "edge": "edge-sublattice channel",
        "cls": "CLS form-factor channel",
        "staggered": "staggered channel",
        "hub": "hub-sublattice channel",
    }
    ticks, tick_labels = tick_positions(manifest, rows)

    figure, axes = plt.subplots(
        len(channels),
        1,
        figsize=(7.2, 1.65 * len(channels) + 0.85),
        sharex=True,
        constrained_layout=True,
    )
    if len(channels) == 1:
        axes = np.asarray([axes])
    for panel, (axis, channel) in enumerate(zip(axes, channels)):
        distance, omega, grid = spectrum_grid(rows, channel, manifest)
        mask = omega <= args.omega_max
        finite = grid[mask][np.isfinite(grid[mask])]
        vmax = float(np.quantile(finite, 0.985)) if finite.size else 1.0
        colormap = mpl.colormaps["magma"].copy()
        colormap.set_bad("#d9d9d9")
        mesh = axis.pcolormesh(
            cell_edges(distance),
            cell_edges(omega[mask]),
            grid[mask],
            shading="flat",
            cmap=colormap,
            vmin=0.0,
            vmax=max(vmax, 1.0e-12),
            rasterized=True,
        )
        lswt_x = np.asarray([float(row["path_distance"]) for row in lswt])
        for branch in range(1, 6):
            lswt_y = np.asarray([float(row[f"omega_sorted_{branch}"]) for row in lswt])
            axis.plot(lswt_x, lswt_y, color="white", linewidth=0.55, alpha=0.82)
        axis.set_ylim(0.0, args.omega_max)
        axis.set_ylabel(r"$\omega/J$")
        axis.set_title(channel_titles.get(channel, channel), loc="left", fontsize=7.5, pad=2)
        axis.text(-0.07, 1.02, chr(ord("a") + panel), transform=axis.transAxes, fontweight="bold", fontsize=8)
        colorbar = figure.colorbar(mesh, ax=axis, pad=0.012, aspect=18)
        colorbar.set_label("SAC intensity", fontsize=6.5)
        colorbar.ax.tick_params(labelsize=6)
        for position in ticks:
            axis.axvline(position, color="white", linewidth=0.35, alpha=0.45)
        missing = np.flatnonzero(np.all(~np.isfinite(grid), axis=0))
        if missing.size:
            missing_label = "elastic\n(omitted)" if channel == "total" else "elastic /\nunresolved"
            for index in missing:
                if index == 0:
                    text_x = distance[index] + 0.015 * (distance[-1] - distance[0])
                    alignment = "left"
                elif index == len(distance) - 1:
                    text_x = distance[index] - 0.015 * (distance[-1] - distance[0])
                    alignment = "right"
                else:
                    text_x = distance[index]
                    alignment = "center"
                axis.text(
                    text_x,
                    0.88 * args.omega_max,
                    missing_label,
                    ha=alignment,
                    va="top",
                    color="#333333",
                    fontsize=5.5,
                )
    axes[-1].set_xticks(ticks, tick_labels)
    axes[-1].set_xlabel(r"momentum path $\Gamma$–M–K–$\Gamma$")
    figure.suptitle(f"Interacting spin spectrum ({args.title_suffix}); white: LSWT", fontsize=8)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix, options in {
        ".svg": {},
        ".pdf": {},
        ".tiff": {"dpi": 600},
        ".png": {"dpi": 240},
    }.items():
        figure.savefig(args.output.with_suffix(suffix), bbox_inches="tight", **options)
    plt.close(figure)


if __name__ == "__main__":
    main()
