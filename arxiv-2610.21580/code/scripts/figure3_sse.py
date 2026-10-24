#!/usr/bin/env python3
"""生成 SSE 收敛与有限尺寸外推图 3。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, NullLocator


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 6.8
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["xtick.major.width"] = 0.7
plt.rcParams["ytick.major.width"] = 0.7
plt.rcParams["xtick.major.size"] = 3.0
plt.rcParams["ytick.major.size"] = 3.0
plt.rcParams["legend.frameon"] = False


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PILOT_CSV = ROOT / "results" / "sse_beta_pilot_v2" / "scan_observables.csv"
DEFAULT_RUN_DIR = ROOT / "results" / "sse_fss_v1"
DEFAULT_OUTPUT = ROOT / "figures" / "figure3_sse"
EXPECTED_SIZES = (4, 6, 8, 10, 12, 14, 16)

DATA_COLOR = "#2F6690"
FIT_COLOR = "#B23A48"
EXACT_COLOR = "#4F5D2F"
SUBLATTICE_A = "#2F6690"
SUBLATTICE_B = "#E09F3E"
PILOT_COLORS = {2: "#B8CBD8", 4: "#7EA8C0", 6: "#477F9F", 8: "#174A67"}
PILOT_MARKERS = {2: "o", 4: "s", 6: "D", 8: "^"}


def fit_files(run_dir: Path) -> dict[str, Path]:
    """生成拟合文件路径。"""
    fit_dir = run_dir / "fits"
    return {
        "energy": fit_dir / "energy_fit.csv",
        "staggered": fit_dir / "staggered_m2_fit.csv",
        "stiffness": fit_dir / "stiffness_fit.csv",
        "plus": fit_dir / "sublattice_plus_m2_secondary_fit.csv",
        "minus": fit_dir / "sublattice_minus_m2_secondary_fit.csv",
    }


def require_inputs(pilot_csv: Path, run_dir: Path, fits: dict[str, Path]) -> None:
    """数据不齐时停止，禁止生成占位图。"""
    required = [
        pilot_csv,
        run_dir / "fss_observables.csv",
        run_dir / "fss_summary.json",
        run_dir / "fit_results_summary.json",
        *fits.values(),
    ]
    missing = [path for path in required if not path.is_file()]
    if missing:
        names = "\n".join(f"  - {path}" for path in missing)
        raise SystemExit(f"Figure 3 数据尚未齐全，未生成图件：\n{names}")


def read_table(path: Path) -> dict[str, np.ndarray]:
    """读取数值表。"""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"空数据表：{path}")
    columns: dict[str, np.ndarray] = {}
    for name in rows[0]:
        try:
            columns[name] = np.asarray([float(row[name]) for row in rows])
        except ValueError:
            continue
    return columns


def read_fit(
    path: Path,
    required: tuple[str, ...],
    pinned_zero: tuple[str, ...],
) -> dict[str, float | tuple[float, float]]:
    """读取固定形式拟合参数。"""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        raw_rows = list(csv.DictReader(handle))
    parameters = [row.get("parameter", "") for row in raw_rows]
    if len(parameters) != len(set(parameters)):
        raise ValueError(f"拟合表存在重复参数：{path}")
    rows = {row["parameter"]: row for row in raw_rows}
    names = (*required, "chi2_per_dof", "failed_resamples")
    missing = [name for name in names if name not in rows]
    if missing:
        raise ValueError(f"拟合表缺少参数 {missing}：{path}")
    try:
        values = {name: float(rows[name]["value"]) for name in names}
        ci = (float(rows["c0"]["ci_lo"]), float(rows["c0"]["ci_hi"]))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"拟合表含非数值字段：{path}") from exc
    if any(not math.isfinite(value) for value in (*values.values(), *ci)):
        raise ValueError(f"拟合表含非有限值：{path}")
    failed = values["failed_resamples"]
    if not failed.is_integer() or failed < 0.0:
        raise ValueError(f"failed_resamples 不是非负整数：{path}")
    if int(failed) != 0:
        raise ValueError(f"bootstrap 存在失败样本：{path}")
    invalid_pins = [name for name in pinned_zero if abs(values[name]) > 1.0e-12]
    if invalid_pins:
        raise ValueError(f"预注册固定系数不是零 {invalid_pins}：{path}")
    if ci[0] > values["c0"] or values["c0"] > ci[1]:
        raise ValueError(f"c0 不在 95% CI 内：{path}")
    if values["chi2_per_dof"] < 0.0:
        raise ValueError(f"chi2_per_dof 不能为负：{path}")
    result: dict[str, float | tuple[float, float]] = {
        name: values[name] for name in required
    }
    result["c0_ci95"] = ci
    result["chi2_per_dof"] = values["chi2_per_dof"]
    return result


def validate_fit_curve(
    fit: dict[str, float | tuple[float, float]],
    fss: dict[str, np.ndarray],
    observable: str,
    coefficients: tuple[str, ...],
    pinned_zero: tuple[str, ...],
) -> None:
    """核对图中曲线与声明的数据列。"""
    prediction = sum(float(fit[name]) / fss["L"] ** int(name[1:]) for name in coefficients)
    residual = (fss[observable] - prediction) / fss[f"{observable}_error"]
    free_count = len(coefficients) - len(pinned_zero)
    chi2 = float(np.sum(residual**2) / max(len(residual) - free_count, 1))
    if not math.isclose(chi2, float(fit["chi2_per_dof"]), rel_tol=1.0e-9, abs_tol=1.0e-10):
        raise ValueError(f"图中拟合与声明数据列不一致：{observable}")


def validate(pilot: dict[str, np.ndarray], fss: dict[str, np.ndarray]) -> None:
    """核对正式尺寸和必需列。"""
    pilot_required = (
        "L",
        "beta_over_L2",
        "beta",
        "energy",
        "energy_error",
        "ferrimagnetic_m2",
        "ferrimagnetic_m2_error",
    )
    fss_required = (
        "L",
        "N",
        "beta",
        "replicas",
        "pooled_bins",
        "energy_per_site",
        "energy_per_site_error",
        "staggered_m2",
        "staggered_m2_error",
        "ferrimagnetic_m2",
        "ferrimagnetic_m2_error",
        "sublattice_plus_m2",
        "sublattice_plus_m2_error",
        "sublattice_minus_m2",
        "sublattice_minus_m2_error",
        "spin_stiffness_cell_average",
        "spin_stiffness_cell_average_error",
    )
    for name in pilot_required:
        if name not in pilot:
            raise ValueError(f"pilot 表缺少列：{name}")
    for name in fss_required:
        if name not in fss:
            raise ValueError(f"FSS 表缺少列：{name}")
    expected_pilot = {(size, ratio) for size in (2, 4, 6, 8) for ratio in (0.25, 0.5, 1.0, 2.0)}
    observed_pilot = list(zip(pilot["L"], np.round(pilot["beta_over_L2"], 8)))
    if len(observed_pilot) != 16 or set(observed_pilot) != expected_pilot:
        raise ValueError("pilot 不是 4×4 beta J/L² 扫描")
    if not np.allclose(pilot["beta"], pilot["beta_over_L2"] * pilot["L"] ** 2, rtol=0.0, atol=1.0e-12):
        raise ValueError("pilot beta 与 beta J/L² 不一致")
    if len(fss["L"]) != 7 or list(fss["L"]) != list(EXPECTED_SIZES):
        raise ValueError("FSS 尺寸轴不完整或未排序")
    for source, names, label in ((pilot, pilot_required, "pilot"), (fss, fss_required, "FSS")):
        for name in names:
            if np.any(~np.isfinite(source[name])):
                raise ValueError(f"{label} 列存在非有限值：{name}")
            if name.endswith("_error") and np.any(source[name] <= 0.0):
                raise ValueError(f"{label} 误差列存在非正数：{name}")
    for size, site_count, beta, replicas, pooled_bins in zip(
        fss["L"], fss["N"], fss["beta"], fss["replicas"], fss["pooled_bins"]
    ):
        if site_count != 5 * size**2 or beta != 2 * size**2 or replicas != 2 or pooled_bins != 128:
            raise ValueError(f"FSS 正式协议字段不一致：L={size:g}")


def validate_summaries(run_dir: Path) -> None:
    """核对汇总元数据。"""
    summary = json.loads((run_dir / "fss_summary.json").read_text(encoding="utf-8"))
    expected = {
        "cell_count": 14,
        "size_count": 7,
        "sizes": list(EXPECTED_SIZES),
        "replicas_per_size": 2,
        "bins_per_size": [128] * 7,
    }
    for name, target in expected.items():
        if summary.get(name) != target:
            raise ValueError(f"FSS 汇总协议字段不一致：{name}")
    fit_summary = json.loads((run_dir / "fit_results_summary.json").read_text(encoding="utf-8"))
    if set(fit_summary.get("primary", {})) != {
        "energy_per_site",
        "staggered_m2",
        "spin_stiffness_cell_average",
    }:
        raise ValueError("拟合汇总主分析集合不完整")


def style_axis(ax: plt.Axes) -> None:
    """统一坐标轴。"""
    ax.tick_params(direction="out", pad=2)
    ax.spines["left"].set_color("#555C62")
    ax.spines["bottom"].set_color("#555C62")


def panel_header(ax: plt.Axes, label: str, title: str) -> None:
    """添加面板标记。"""
    ax.text(-0.17, 1.025, label, transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="bottom")
    ax.set_title(title, loc="left", fontsize=7.8, pad=4)


def error_points(ax: plt.Axes, x: np.ndarray, y: np.ndarray, error: np.ndarray, **kwargs) -> None:
    """绘制统一误差棒。"""
    defaults = {
        "color": DATA_COLOR,
        "marker": "o",
        "markersize": 3.6,
        "markerfacecolor": "white",
        "markeredgewidth": 0.85,
        "linestyle": "none",
        "elinewidth": 0.75,
        "capsize": 2.0,
        "capthick": 0.75,
        "zorder": 3,
    }
    defaults.update(kwargs)
    ax.errorbar(x, y, yerr=error, **defaults)


def fit_note(ax: plt.Axes, fit: dict[str, float | tuple[float, float]], symbol: str) -> None:
    """标注外推值和拟合质量。"""
    c0 = float(fit["c0"])
    ci_lo, ci_hi = fit["c0_ci95"]
    chi2 = float(fit["chi2_per_dof"])
    text = rf"$({symbol})_\infty={c0:.5g}$" + "\n" + rf"95% CI [{ci_lo:.5g}, {ci_hi:.5g}]" + "\n" + rf"$\chi^2_\nu={chi2:.2f}$"
    ax.text(
        0.035,
        0.045,
        text,
        transform=ax.transAxes,
        fontsize=5.5,
        va="bottom",
        color="#3D454B",
        bbox={"facecolor": "white", "edgecolor": "#D1D6DA", "linewidth": 0.5, "alpha": 0.90, "pad": 1.3},
        zorder=5,
    )


def inverse_square_axis(ax: plt.Axes) -> None:
    """设置 1/L² 横轴。"""
    ax.set_xlim(-0.002, 0.065)
    ax.set_xticks([0.0, 0.02, 0.04, 0.06])
    ax.set_xlabel(r"$1/L^2$", labelpad=2)
    style_axis(ax)


def draw_pilot(panel, pilot: dict[str, np.ndarray]) -> None:
    """绘制 beta J/L² 收敛。"""
    subgrid = panel.subgridspec(2, 1, hspace=0.08)
    ax_e = plt.gcf().add_subplot(subgrid[0])
    ax_m = plt.gcf().add_subplot(subgrid[1], sharex=ax_e)
    for size in (2, 4, 6, 8):
        select = pilot["L"].astype(int) == size
        order = np.argsort(pilot["beta_over_L2"][select])
        x = pilot["beta_over_L2"][select][order]
        common = {
            "color": PILOT_COLORS[size],
            "marker": PILOT_MARKERS[size],
            "markersize": 2.8,
            "linewidth": 0.9,
            "elinewidth": 0.55,
            "capsize": 1.4,
        }
        ax_e.errorbar(x, pilot["energy"][select][order], yerr=pilot["energy_error"][select][order], label=rf"$L={size}$", **common)
        ax_m.errorbar(
            x,
            pilot["ferrimagnetic_m2"][select][order],
            yerr=pilot["ferrimagnetic_m2_error"][select][order],
            **common,
        )
    for ax in (ax_e, ax_m):
        ax.set_xscale("log", base=2)
        ax.set_xlim(0.22, 2.25)
        ax.set_xticks([0.25, 0.5, 1.0, 2.0])
        ax.xaxis.set_minor_locator(NullLocator())
        style_axis(ax)
    ax_e.tick_params(labelbottom=False)
    ax_e.set_ylabel(r"$E/(NJ)$", labelpad=2)
    ax_m.set_ylabel(r"$m_f^2$", labelpad=2)
    ax_m.set_xlabel(r"$\beta J/L^2$", labelpad=1)
    ax_m.set_xticklabels(["0.25", "0.5", "1", "2"])
    ax_e.legend(loc="lower right", ncol=2, fontsize=5.3, columnspacing=0.8, handlelength=1.5, handletextpad=0.35)
    panel_header(ax_e, "a", r"Fixed-aspect-ratio pilot")


def draw_energy(ax: plt.Axes, fss: dict[str, np.ndarray], fit: dict[str, float | tuple[float, float]]) -> None:
    """绘制能量主拟合。"""
    sizes = fss["L"]
    x = 1.0e3 / sizes**4
    error_points(ax, x, fss["energy_per_site"], fss["energy_per_site_error"])
    grid = np.linspace(0.0, max(x) * 1.025, 300)
    u = grid / 1.0e3
    curve = float(fit["c0"]) + float(fit["c4"]) * u + float(fit["c6"]) * u**1.5
    ax.plot(grid, curve, color=FIT_COLOR, linewidth=1.35, label="primary fit")
    ax.set_xlim(-0.10, max(x) * 1.06)
    ax.set_xlabel(r"$10^3/L^4$", labelpad=2)
    ax.set_ylabel(r"Energy per site  $E/(NJ)$")
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.5f"))
    panel_header(ax, "b", r"Energy: $L^{-4}+L^{-6}$")
    fit_note(ax, fit, "E/(NJ)")
    style_axis(ax)


def draw_even_fit(
    ax: plt.Axes,
    fss: dict[str, np.ndarray],
    fit: dict[str, float | tuple[float, float]],
    value_name: str,
    error_name: str,
    ylabel: str,
    label: str,
    title: str,
    symbol: str,
) -> None:
    """绘制偶次 1/L 外推。"""
    x = 1.0 / fss["L"] ** 2
    error_points(ax, x, fss[value_name], fss[error_name])
    grid = np.linspace(0.0, max(x) * 1.025, 300)
    curve = float(fit["c0"]) + float(fit["c2"]) * grid + float(fit["c4"]) * grid**2
    ax.plot(grid, curve, color=FIT_COLOR, linewidth=1.35)
    ax.set_ylabel(ylabel)
    inverse_square_axis(ax)
    panel_header(ax, label, title)
    fit_note(ax, fit, symbol)


def draw_ferrimagnetic(ax: plt.Axes, fss: dict[str, np.ndarray]) -> None:
    """绘制 Lieb--Mattis 精确约束。"""
    x = 1.0 / fss["L"] ** 2
    error_points(ax, x, fss["ferrimagnetic_m2"], fss["ferrimagnetic_m2_error"])
    grid = np.linspace(0.0, max(x) * 1.025, 300)
    ax.plot(grid, 0.01 + 0.02 * grid, color=EXACT_COLOR, linewidth=1.35, linestyle=(0, (4, 2)), label="exact")
    ax.legend(loc="best", fontsize=5.8, handlelength=2.0)
    ax.set_ylabel(r"Ferrimagnetic moment  $m_f^2$")
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    inverse_square_axis(ax)
    panel_header(ax, "d", "Exact net-moment constraint")
    ax.text(0.035, 0.06, r"$m_f^2=0.01+0.02/L^2$", transform=ax.transAxes, fontsize=5.8, color=EXACT_COLOR)


def draw_sublattices(
    ax: plt.Axes,
    fss: dict[str, np.ndarray],
    fit_a: dict[str, float | tuple[float, float]],
    fit_b: dict[str, float | tuple[float, float]],
) -> None:
    """绘制两个子晶格矩。"""
    x = 1.0 / fss["L"] ** 2
    error_points(
        ax,
        x,
        fss["sublattice_plus_m2"],
        fss["sublattice_plus_m2_error"],
        color=SUBLATTICE_A,
        marker="o",
        label=r"$m_A^2$",
    )
    error_points(
        ax,
        x,
        fss["sublattice_minus_m2"],
        fss["sublattice_minus_m2_error"],
        color=SUBLATTICE_B,
        marker="s",
        label=r"$m_B^2$",
    )
    grid = np.linspace(0.0, max(x) * 1.025, 300)
    for fit, color in ((fit_a, SUBLATTICE_A), (fit_b, SUBLATTICE_B)):
        curve = float(fit["c0"]) + float(fit["c2"]) * grid + float(fit["c4"]) * grid**2
        ax.plot(grid, curve, color=color, linewidth=1.25)
    ax.legend(loc="best", fontsize=5.8, ncol=2, columnspacing=0.9, handletextpad=0.35)
    ax.set_ylabel("Squared sublattice moment")
    inverse_square_axis(ax)
    panel_header(ax, "e", "Sublattice moments")
    ax.text(0.035, 0.055, "secondary fits", transform=ax.transAxes, fontsize=5.6, color="#596269")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--pilot-csv", type=Path, default=DEFAULT_PILOT_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    run_dir = args.run_dir.resolve()
    pilot_csv = args.pilot_csv.resolve()
    output = args.output.resolve()
    files = fit_files(run_dir)
    require_inputs(pilot_csv, run_dir, files)
    validate_summaries(run_dir)
    pilot = read_table(pilot_csv)
    fss = read_table(run_dir / "fss_observables.csv")
    validate(pilot, fss)
    fits = {
        "energy": read_fit(
            files["energy"],
            ("c0", "c1", "c2", "c3", "c4", "c5", "c6"),
            ("c1", "c2", "c3", "c5"),
        ),
        "staggered": read_fit(
            files["staggered"], ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")
        ),
        "stiffness": read_fit(
            files["stiffness"], ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")
        ),
        "plus": read_fit(
            files["plus"], ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")
        ),
        "minus": read_fit(
            files["minus"], ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")
        ),
    }
    fit_specs = {
        "energy": ("energy_per_site", ("c0", "c1", "c2", "c3", "c4", "c5", "c6"), ("c1", "c2", "c3", "c5")),
        "staggered": ("staggered_m2", ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")),
        "stiffness": ("spin_stiffness_cell_average", ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")),
        "plus": ("sublattice_plus_m2", ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")),
        "minus": ("sublattice_minus_m2", ("c0", "c1", "c2", "c3", "c4"), ("c1", "c3")),
    }
    for name, (observable, coefficients, pins) in fit_specs.items():
        validate_fit_curve(fits[name], fss, observable, coefficients, pins)

    width_mm, height_mm = 183.0, 132.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    grid = fig.add_gridspec(
        2,
        3,
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.965,
        wspace=0.38,
        hspace=0.40,
    )
    draw_pilot(grid[0, 0], pilot)
    draw_energy(fig.add_subplot(grid[0, 1]), fss, fits["energy"])
    draw_even_fit(
        fig.add_subplot(grid[0, 2]),
        fss,
        fits["staggered"],
        "staggered_m2",
        "staggered_m2_error",
        r"Staggered moment  $m_s^2$",
        "c",
        "Antiferromagnetic order",
        "m_s^2",
    )
    draw_ferrimagnetic(fig.add_subplot(grid[1, 0]), fss)
    draw_sublattices(fig.add_subplot(grid[1, 1]), fss, fits["plus"], fits["minus"])
    draw_even_fit(
        fig.add_subplot(grid[1, 2]),
        fss,
        fits["stiffness"],
        "spin_stiffness_cell_average",
        "spin_stiffness_cell_average_error",
        r"Cell-axis stiffness  $\rho_{\rm cell}/J$",
        "f",
        "Spin stiffness",
        r"\rho_{\rm cell}/J",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".svg"), facecolor="white")
    fig.savefig(output.with_suffix(".pdf"), facecolor="white")
    fig.savefig(output.with_suffix(".tiff"), dpi=600, facecolor="white")
    fig.savefig(output.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Figure 3 已写入：{output.parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
