#!/usr/bin/env python3
"""校验正式 SSE 数据并安全回填论文。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

try:
    from scripts.analyze_sse_fss import summarize as recompute_fss
except ModuleNotFoundError:  # 直接运行脚本时使用
    from analyze_sse_fss import summarize as recompute_fss


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SIZES = (4, 6, 8, 10, 12, 14, 16)
EXPECTED_FIGURE_SUFFIXES = (".pdf", ".svg", ".png", ".tiff")
MARKER = "SSE_PRODUCTION_BACKFILL_V1"

MD_ABSTRACT_PLACEHOLDER = (
    "[TODO: insert the final thermodynamic-limit energy, staggered moment, "
    "and cell-axis spin stiffness]"
)
MD_BODY_PLACEHOLDER = (
    "> **Production SSE placeholder.** [TODO: insert pooled observables, "
    "autocorrelation diagnostics, committed finite-size fits, confidence intervals, "
    "and a figure after all 14 sse_fss_v1 cells pass provenance and completion checks. "
    "Do not replace this block with single-replica or incomplete-run values.]"
)
MD_FIGURE_PLACEHOLDER = (
    "**Figure 3 placeholder.** (a) Pilot convergence versus $\\beta J/L^2$. "
    "(b)--(f) Two-replica finite-size data and committed extrapolations for $E/(NJ)$, "
    "$m_s^2$, $m_f^2$, sublattice moments, and $\\rho_{\\rm cell}$. "
    "[TODO: replace after the complete production analysis.]"
)
MD_CONCLUSION_PLACEHOLDER = "[TODO: add the final SSE thermodynamic-limit estimates]"

TEX_ABSTRACT_PLACEHOLDER = (
    r"\todo{insert the final thermodynamic-limit energy, staggered moment, "
    r"and cell-axis spin stiffness}"
)
TEX_BODY_PLACEHOLDER = "\n".join(
    (
        r"\begin{quote}",
        r"\todo{Production SSE placeholder. Insert pooled observables, autocorrelation "
        r"diagnostics, committed finite-size fits, confidence intervals, and a figure "
        r"after all 14 \texttt{sse\_fss\_v1} cells pass provenance and completion checks. "
        r"Do not replace this block with single-replica or incomplete-run values.}",
        r"\end{quote}",
    )
)
TEX_FIGURE_PLACEHOLDER = "\n".join(
    (
        r"\begin{figure*}[t]",
        r"\centering",
        r"\fbox{\parbox{0.94\textwidth}{\centering\vspace{18mm}",
        r"\todo{Figure 3 is pending the complete production SSE analysis.}\\[1ex]",
        r"(a) Pilot convergence versus $\beta J/L^2$. (b)--(f) Two-replica finite-size "
        r"data and committed extrapolations for $E/(NJ)$, $m_s^2$, $m_f^2$, sublattice "
        r"moments, and $\rho_{\rm cell}$.",
        r"\vspace{18mm}}}",
        r"\caption{SSE convergence and finite-size scaling. \todo{Replace this placeholder "
        r"after the complete production analysis.}}",
        r"\label{fig:sse}",
        r"\end{figure*}",
    )
)
TEX_CONCLUSION_PLACEHOLDER = r"\todo{add the final SSE thermodynamic-limit estimates}"

FSS_REQUIRED_COLUMNS = (
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
    "lieb_mattis_target",
    "sublattice_plus_m2",
    "sublattice_plus_m2_error",
    "sublattice_minus_m2",
    "sublattice_minus_m2_error",
    "sublattice_cross",
    "sublattice_cross_error",
    "spin_stiffness_cell_average",
    "spin_stiffness_cell_average_error",
)
PRIMARY_NAMES = (
    "energy_per_site",
    "staggered_m2",
    "spin_stiffness_cell_average",
)
SECONDARY_NAMES = (
    "sublattice_plus_m2",
    "sublattice_minus_m2",
    "sublattice_cross",
)


class BackfillError(RuntimeError):
    """回填门禁错误。"""


def finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise BackfillError(f"{label} 不是数值")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise BackfillError(f"{label} 不是数值") from exc
    if not math.isfinite(result):
        raise BackfillError(f"{label} 不是有限值")
    return result


def reject_nonfinite(value: Any, label: str) -> None:
    """递归拒绝 JSON 非有限数。"""
    if isinstance(value, float) and not math.isfinite(value):
        raise BackfillError(f"{label} 含非有限值")
    if isinstance(value, Mapping):
        for key, item in value.items():
            reject_nonfinite(item, f"{label}.{key}")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for index, item in enumerate(value):
            reject_nonfinite(item, f"{label}[{index}]")


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise BackfillError(f"缺少文件：{path}")
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)),
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise BackfillError(f"JSON 无效或含非有限数：{path}") from exc
    if not isinstance(value, dict):
        raise BackfillError(f"JSON 顶层必须是对象：{path}")
    reject_nonfinite(value, path.name)
    return value


def load_fss_csv(path: Path) -> list[dict[str, float]]:
    if not path.is_file():
        raise BackfillError(f"缺少文件：{path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    if len(fields) != len(set(fields)):
        raise BackfillError("FSS 表存在重复列名")
    missing = [name for name in FSS_REQUIRED_COLUMNS if name not in fields]
    if missing:
        raise BackfillError(f"FSS 表缺少列：{missing}")
    if len(rows) != len(EXPECTED_SIZES):
        raise BackfillError("FSS 表必须恰含七个尺寸")
    numeric: list[dict[str, float]] = []
    for row_index, row in enumerate(rows, start=2):
        if None in row:
            raise BackfillError(f"FSS 第 {row_index} 行字段数量不一致")
        converted = {
            name: finite_float(value, f"FSS 第 {row_index} 行 {name}")
            for name, value in row.items()
        }
        numeric.append(converted)
    validate_fss_rows(numeric)
    return numeric


def validate_fss_rows(rows: list[dict[str, float]]) -> None:
    sizes = [int(row["L"]) for row in rows]
    if any(row["L"] != int(row["L"]) for row in rows) or sizes != list(EXPECTED_SIZES):
        raise BackfillError("FSS 尺寸轴不完整或未排序")
    for row, size in zip(rows, EXPECTED_SIZES):
        expected = {
            "N": 5 * size * size,
            "beta": 2 * size * size,
            "replicas": 2,
            "pooled_bins": 128,
        }
        for name, target in expected.items():
            if row[name] != target:
                raise BackfillError(f"FSS 协议字段不一致：L={size}, {name}")
        target = 0.01 + 0.02 / size**2
        if not math.isclose(row["lieb_mattis_target"], target, rel_tol=0.0, abs_tol=1.0e-12):
            raise BackfillError(f"Lieb--Mattis 目标不一致：L={size}")
        for name, value in row.items():
            if name.endswith("_error") and value <= 0.0:
                raise BackfillError(f"FSS 误差必须为正：L={size}, {name}")


def validate_fss_summary(summary: dict[str, Any]) -> None:
    expected = {
        "cell_count": 14,
        "size_count": 7,
        "sizes": list(EXPECTED_SIZES),
        "replicas_per_size": 2,
        "bins_per_size": [128] * 7,
    }
    for name, target in expected.items():
        if summary.get(name) != target:
            raise BackfillError(f"FSS 汇总协议字段不一致：{name}")
    for name in (
        "max_abs_replica_lag1",
        "max_abs_recomputed_replica_lag1",
        "max_replica_difference_sigma",
        "max_blocking_error_ratio",
    ):
        value = finite_float(summary.get(name), f"fss_summary.{name}")
        if value < 0.0:
            raise BackfillError(f"FSS 诊断量不能为负：{name}")


def fit_node(container: Mapping[str, Any], key: str, label: str) -> dict[str, Any]:
    value = container.get(key)
    if not isinstance(value, dict):
        raise BackfillError(f"缺少拟合结果：{label}")
    center = finite_float(value.get("intercept"), f"{label}.intercept")
    ci = value.get("ci95")
    if not isinstance(ci, list) or len(ci) != 2:
        raise BackfillError(f"{label}.ci95 必须含两个端点")
    lower, upper = (finite_float(ci[0], f"{label}.ci95[0]"), finite_float(ci[1], f"{label}.ci95[1]"))
    if not lower < upper or not lower <= center <= upper:
        raise BackfillError(f"{label} 的中心值与 95% CI 不一致")
    chi2 = finite_float(value.get("chi2_per_dof"), f"{label}.chi2_per_dof")
    if chi2 < 0.0:
        raise BackfillError(f"{label}.chi2_per_dof 不能为负")
    failed = value.get("failed_resamples")
    if isinstance(failed, bool) or not isinstance(failed, int) or failed != 0:
        raise BackfillError(f"{label} 的 bootstrap 存在失败或字段无效")
    return value


def validate_interval(value: Any, center: float, label: str) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise BackfillError(f"{label} 必须含两个端点")
    lower, upper = finite_float(value[0], f"{label}[0]"), finite_float(value[1], f"{label}[1]")
    if not lower < upper or not lower <= center <= upper:
        raise BackfillError(f"{label} 与中心值不一致")
    return lower, upper


def validate_fit_summary(
    result: dict[str, Any], rows: list[dict[str, float]], fss_summary: dict[str, Any]
) -> None:
    primary = result.get("primary")
    secondary = result.get("secondary")
    sensitivity = result.get("sensitivity")
    if not all(isinstance(item, dict) for item in (primary, secondary, sensitivity)):
        raise BackfillError("拟合汇总缺少主拟合、补充拟合或敏感性分析")
    primary_fits = {
        name: fit_node(primary, name, f"primary.{name}") for name in PRIMARY_NAMES
    }
    for name in SECONDARY_NAMES:
        fit_node(secondary, name, f"secondary.{name}")
    for name in PRIMARY_NAMES:
        group = sensitivity.get(name)
        if not isinstance(group, dict):
            raise BackfillError(f"缺少敏感性结果：{name}")
        fits = [primary_fits[name]]
        for minimum in (6, 8):
            key = f"Lmin_{minimum}"
            fits.append(fit_node(group, key, f"sensitivity.{name}.{key}"))
        centers = [finite_float(item["intercept"], name) for item in fits]
        span = group.get("intercept_span")
        if not isinstance(span, list) or len(span) != 2:
            raise BackfillError(f"敏感性范围无效：{name}")
        endpoints = [finite_float(item, f"{name}.intercept_span") for item in span]
        if not all(
            math.isclose(actual, target, rel_tol=1.0e-12, abs_tol=1.0e-14)
            for actual, target in zip(endpoints, (min(centers), max(centers)))
        ):
            raise BackfillError(f"敏感性范围与拟合截距不一致：{name}")

    constrained = result.get(
        "constrained_sublattice_assuming_collinear_order_and_cross_factorization"
    )
    if not isinstance(constrained, dict):
        raise BackfillError("缺少约束子晶格变换")
    constrained_names = (
        "staggered_moment",
        "plus_moment",
        "minus_moment",
        "plus_m2",
        "minus_m2",
        "cross",
    )
    constrained_ci = constrained.get("ci95")
    if not isinstance(constrained_ci, dict):
        raise BackfillError("缺少约束子晶格区间")
    for name in constrained_names:
        center = finite_float(constrained.get(name), f"constrained.{name}")
        validate_interval(constrained_ci.get(name), center, f"constrained.ci95.{name}")
    ms2 = finite_float(primary_fits["staggered_m2"]["intercept"], "primary.staggered_m2")
    staggered_moment = finite_float(
        constrained["staggered_moment"], "constrained.staggered_moment"
    )
    if not math.isclose(staggered_moment**2, ms2, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise BackfillError("约束交错矩与主拟合不一致")

    exact = result.get("ferrimagnetic_exact_check")
    if not isinstance(exact, dict) or exact.get("target") != "0.01+0.02/L^2":
        raise BackfillError("Lieb--Mattis 检验定义不一致")
    residual = max(
        abs((row["ferrimagnetic_m2"] - row["lieb_mattis_target"]) / row["ferrimagnetic_m2_error"])
        for row in rows
    )
    recorded = finite_float(
        exact.get("max_abs_standardized_residual"),
        "ferrimagnetic_exact_check.max_abs_standardized_residual",
    )
    if not math.isclose(recorded, residual, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise BackfillError("Lieb--Mattis 最大标准化残差与 FSS 表不一致")

    diagnostics = result.get("bin_diagnostics")
    if not isinstance(diagnostics, dict):
        raise BackfillError("缺少 bin 诊断汇总")
    mapping = {
        "max_abs_replica_lag1": "max_abs_replica_lag1",
        "max_blocking_error_ratio": "max_blocking_error_ratio",
        "max_replica_difference_sigma": "max_replica_difference_sigma",
    }
    for fit_name, fss_name in mapping.items():
        actual = finite_float(diagnostics.get(fit_name), f"bin_diagnostics.{fit_name}")
        target = finite_float(fss_summary.get(fss_name), f"fss_summary.{fss_name}")
        if not math.isclose(actual, target, rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise BackfillError(f"拟合与 FSS 诊断量不一致：{fit_name}")

    poor = [
        name
        for name, fit in primary_fits.items()
        if finite_float(fit["chi2_per_dof"], name) > 3.0
    ]
    if poor:
        raise BackfillError(f"主拟合 chi2/dof>3，需人工决定降级措辞：{poor}")


def compare_recomputed(
    rows: list[dict[str, float]],
    summary: dict[str, Any],
    recomputed_rows: list[dict[str, Any]],
    recomputed_summary: dict[str, Any],
) -> None:
    """以原始清单和 bin 重算结果复核汇总。"""
    if len(recomputed_rows) != len(rows):
        raise BackfillError("原始单元重算后的尺寸数不一致")
    for row, fresh in zip(rows, recomputed_rows):
        for name, value in row.items():
            if name not in fresh:
                raise BackfillError(f"原始单元重算缺少列：{name}")
            target = finite_float(fresh[name], f"重算 {name}")
            if not math.isclose(value, target, rel_tol=1.0e-12, abs_tol=1.0e-14):
                raise BackfillError(f"FSS 表与原始单元重算不一致：L={row['L']:g}, {name}")
    for name in (
        "cell_count",
        "size_count",
        "sizes",
        "replicas_per_size",
        "bins_per_size",
    ):
        if summary.get(name) != recomputed_summary.get(name):
            raise BackfillError(f"FSS 汇总与原始单元重算不一致：{name}")
    for name in (
        "max_abs_replica_lag1",
        "max_abs_recomputed_replica_lag1",
        "max_replica_difference_sigma",
        "max_blocking_error_ratio",
    ):
        actual = finite_float(summary.get(name), f"fss_summary.{name}")
        target = finite_float(recomputed_summary.get(name), f"重算 {name}")
        if not math.isclose(actual, target, rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise BackfillError(f"FSS 诊断量与原始单元重算不一致：{name}")


def require_figures(figure_base: Path) -> None:
    missing = [figure_base.with_suffix(suffix) for suffix in EXPECTED_FIGURE_SUFFIXES if not figure_base.with_suffix(suffix).is_file()]
    if missing:
        raise BackfillError("Figure 3 图件不完整：" + ", ".join(str(path) for path in missing))
    if any(path.stat().st_size == 0 for path in (figure_base.with_suffix(s) for s in EXPECTED_FIGURE_SUFFIXES)):
        raise BackfillError("Figure 3 存在空文件")


def validate_placeholders(markdown: str, tex: str) -> None:
    if MARKER in markdown or MARKER in tex:
        raise BackfillError("论文已经完成正式 SSE 回填，拒绝二次应用")
    for label, text, placeholder in (
        ("Markdown 摘要", markdown, MD_ABSTRACT_PLACEHOLDER),
        ("Markdown 正文", markdown, MD_BODY_PLACEHOLDER),
        ("Markdown Figure 3", markdown, MD_FIGURE_PLACEHOLDER),
        ("Markdown 结论", markdown, MD_CONCLUSION_PLACEHOLDER),
        ("LaTeX 摘要", tex, TEX_ABSTRACT_PLACEHOLDER),
        ("LaTeX 正文", tex, TEX_BODY_PLACEHOLDER),
        ("LaTeX Figure 3", tex, TEX_FIGURE_PLACEHOLDER),
        ("LaTeX 结论", tex, TEX_CONCLUSION_PLACEHOLDER),
    ):
        count = text.count(placeholder)
        if count != 1:
            raise BackfillError(f"{label} 占位符数量应为 1，实际为 {count}")


def fit_numbers(fit: Mapping[str, Any]) -> tuple[str, str, str]:
    center = finite_float(fit["intercept"], "intercept")
    lower, upper = (finite_float(fit["ci95"][0], "ci95[0]"), finite_float(fit["ci95"][1], "ci95[1]"))
    half_width = max(center - lower, upper - center)
    decimals = max(0, 1 - math.floor(math.log10(half_width)))
    return tuple(f"{value:.{decimals}f}" for value in (center, lower, upper))


def plain_number(value: Any, digits: int = 3) -> str:
    return f"{finite_float(value, 'value'):.{digits}g}"


def build_replacements(result: dict[str, Any]) -> tuple[dict[str, str], dict[str, str]]:
    primary = result["primary"]
    secondary = result["secondary"]
    sensitivity = result["sensitivity"]
    constrained = result[
        "constrained_sublattice_assuming_collinear_order_and_cross_factorization"
    ]
    diagnostics = result["bin_diagnostics"]

    e0, e0_lo, e0_hi = fit_numbers(primary["energy_per_site"])
    ms2, ms2_lo, ms2_hi = fit_numbers(primary["staggered_m2"])
    rho, rho_lo, rho_hi = fit_numbers(primary["spin_stiffness_cell_average"])
    ms, ms_lo, ms_hi = fit_numbers(
        {
            "intercept": constrained["staggered_moment"],
            "ci95": constrained["ci95"]["staggered_moment"],
        }
    )
    m_a, m_a_lo, m_a_hi = fit_numbers(
        {"intercept": constrained["plus_moment"], "ci95": constrained["ci95"]["plus_moment"]}
    )
    m_b, m_b_lo, m_b_hi = fit_numbers(
        {"intercept": constrained["minus_moment"], "ci95": constrained["ci95"]["minus_moment"]}
    )
    secondary_text: dict[str, tuple[str, str, str, str]] = {}
    for name in SECONDARY_NAMES:
        center, lower, upper = fit_numbers(secondary[name])
        secondary_text[name] = (
            center,
            lower,
            upper,
            f"{finite_float(secondary[name]['chi2_per_dof'], name):.2f}",
        )

    lag1 = plain_number(diagnostics["max_abs_replica_lag1"])
    replica_z = plain_number(diagnostics["max_replica_difference_sigma"])
    blocking = plain_number(diagnostics["max_blocking_error_ratio"])
    mf_residual = plain_number(
        result["ferrimagnetic_exact_check"]["max_abs_standardized_residual"]
    )
    e0_span = "[" + ", ".join(
        plain_number(value, 6) for value in sensitivity["energy_per_site"]["intercept_span"]
    ) + "]"
    ms2_span = "[" + ", ".join(
        plain_number(value, 6) for value in sensitivity["staggered_m2"]["intercept_span"]
    ) + "]"
    rho_span = "[" + ", ".join(
        plain_number(value, 6)
        for value in sensitivity["spin_stiffness_cell_average"]["intercept_span"]
    ) + "]"
    chi2 = {
        name: f"{finite_float(primary[name]['chi2_per_dof'], name):.2f}"
        for name in PRIMARY_NAMES
    }
    a2, b2, cross = (secondary_text[name] for name in SECONDARY_NAMES)

    abstract_md = (
        f"fixed-aspect-ratio thermodynamic-limit estimates $E_0/(NJ)={e0}$ "
        f"[95% bootstrap CI ${e0_lo}$--${e0_hi}$], $m_{{s,\\infty}}={ms}$ "
        f"[${ms_lo}$--${ms_hi}$], and $\\rho_{{{{\\rm cell}},\\infty}}/J={rho}$ "
        f"[${rho_lo}$--${rho_hi}$]"
    )
    abstract_tex = abstract_md.replace("95%", r"95\%")

    body_md = f"""<!-- {MARKER} -->
All 14 production cells passed the completion, provenance, and artifact-hash checks. Each size therefore contributes 128 pooled bins from two independent replicas. Across the production set, the largest absolute bin lag-one autocorrelation was ${lag1}$, the largest replica difference was ${replica_z}\\,\\sigma$, and reblocking increased an individual standard error by at most a factor of ${blocking}$.

The preregistered $L=4$--16 fits gave $E_0/(NJ)={e0}$ with 95% bootstrap CI $[{e0_lo},{e0_hi}]$ and $\\chi^2/{{\\rm dof}}={chi2['energy_per_site']}$, $m_{{s,\\infty}}^2={ms2}$ with $[{ms2_lo},{ms2_hi}]$ and $\\chi^2/{{\\rm dof}}={chi2['staggered_m2']}$, and $\\rho_{{{{\\rm cell}},\\infty}}/J={rho}$ with $[{rho_lo},{rho_hi}]$ and $\\chi^2/{{\\rm dof}}={chi2['spin_stiffness_cell_average']}$. Applying the same forms at $L_{{\\min}}=6$ and 8 gave intercept ranges ${e0_span}$, ${ms2_span}$, and ${rho_span}$, respectively. These ranges measure fit-window sensitivity; they are not statistical intervals.

The finite-size ferrimagnetic estimator was compared with the zero-temperature Lieb--Mattis target $m_f^2=0.01+0.02/L^2$; the maximum standardized residual was ${mf_residual}$. Direct secondary fits gave $m_A^2={a2[0]}$ [95% CI ${a2[1]}$--${a2[2]}$, $\\chi^2/{{\\rm dof}}={a2[3]}$], $m_B^2={b2[0]}$ [${b2[1]}$--${b2[2]}$, $\\chi^2/{{\\rm dof}}={b2[3]}$], and $\\langle\\mathbf M_A\\!\\cdot\\!\\mathbf M_B\\rangle/(N_A N_B)={cross[0]}$ [${cross[1]}$--${cross[2]}$, $\\chi^2/{{\\rm dof}}={cross[3]}$]. Under the additional assumptions of collinear antiparallel order and thermodynamic cross-term factorization, the primary staggered fit instead implies $m_A={m_a}$ [${m_a_lo}$--${m_a_hi}$] and $m_B={m_b}$ [${m_b_lo}$--${m_b_hi}$]. We report this constrained conversion separately from the direct secondary fits.

Error bars are the largest of the unblocked and replica-aware 2/4/8-bin standard errors. The quoted confidence intervals are statistical intervals for the registered ansatzes; they do not include systematic uncertainty from the fit form or the finite fixed aspect ratio $\\beta J/L^2=2$."""

    caption_md = (
        "**Figure 3. SSE convergence and fixed-aspect-ratio finite-size scaling.** "
        "(a) Pilot energy and squared ferrimagnetic moment versus $\\beta J/L^2$ for "
        "$L=2,4,6,8$. (b) Energy per site fitted to $e_0+a_4L^{-4}+a_6L^{-6}$. "
        "(c) Squared staggered moment fitted to $m_{s,\\infty}^2+b_2L^{-2}+b_4L^{-4}$. "
        "(d) Squared ferrimagnetic moment compared with the finite-size zero-temperature "
        "Lieb--Mattis target $0.01+0.02/L^2$. (e) Direct secondary fits of the hub and "
        "edge squared sublattice moments. (f) Cell-axis spin stiffness fitted to "
        "$\\rho_{{\\rm cell},\\infty}+c_2L^{-2}+c_4L^{-4}$. Production points in "
        "(b)--(f) pool two replicas and 128 bins per size; error bars are the largest of "
        "the unblocked and replica-aware 2/4/8-bin standard errors. Curves in (b), (c), "
        "and (f) are the preregistered primary fits, whereas those in (e) are secondary. "
        "In-panel intervals are 95% parametric-bootstrap statistical intervals for the "
        "stated fixed forms."
    )
    figure_md = f"![SSE convergence and fixed-aspect-ratio finite-size scaling.](../figures/figure3_sse.png)\n\n{caption_md}"
    conclusion_md = (
        f"fixed-aspect-ratio SSE estimates $E_0/(NJ)={e0}$, "
        f"$m_{{s,\\infty}}={ms}$, and $\\rho_{{{{\\rm cell}},\\infty}}/J={rho}$, "
        "with the finite-size net moment tested against the Lieb--Mattis constraint"
    )

    body_tex = body_md.replace("95%", r"95\%").replace(
        f"<!-- {MARKER} -->", f"% {MARKER}"
    )
    caption_tex = caption_md.replace(
        "**Figure 3. SSE convergence and fixed-aspect-ratio finite-size scaling.** ", ""
    ).replace("95%", r"95\%")
    figure_tex = "\n".join(
        (
            r"\begin{figure*}[t]",
            r"\centering",
            r"\includegraphics[width=\textwidth]{../figures/figure3_sse.pdf}",
            rf"\caption{{{caption_tex}}}",
            r"\label{fig:sse}",
            r"\end{figure*}",
        )
    )
    return (
        {
            MD_ABSTRACT_PLACEHOLDER: abstract_md,
            MD_BODY_PLACEHOLDER: body_md,
            MD_FIGURE_PLACEHOLDER: figure_md,
            MD_CONCLUSION_PLACEHOLDER: conclusion_md,
        },
        {
            TEX_ABSTRACT_PLACEHOLDER: abstract_tex,
            TEX_BODY_PLACEHOLDER: body_tex,
            TEX_FIGURE_PLACEHOLDER: figure_tex,
            TEX_CONCLUSION_PLACEHOLDER: conclusion_md,
        },
    )


def replace_once(text: str, replacements: Mapping[str, str]) -> str:
    for old, new in replacements.items():
        if text.count(old) != 1:
            raise BackfillError("占位符在写入前发生变化")
        text = text.replace(old, new)
    return text


def backup_path(path: Path) -> Path:
    return path.with_name(path.name + ".pre_sse_backfill.bak")


def stage_file(path: Path, content: bytes) -> Path:
    descriptor, raw_path = tempfile.mkstemp(prefix=f".{path.name}.sse-", dir=path.parent)
    staged = Path(raw_path)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        staged.unlink(missing_ok=True)
        raise
    return staged


def apply_transaction(files: Mapping[Path, str]) -> list[Path]:
    originals = {path: path.read_bytes() for path in files}
    backups = {path: backup_path(path) for path in files}
    existing = [path for path in backups.values() if path.exists()]
    if existing:
        raise BackfillError("备份文件已存在，拒绝覆盖：" + ", ".join(str(path) for path in existing))
    staged = {path: stage_file(path, content.encode("utf-8")) for path, content in files.items()}
    created_backups: list[Path] = []
    replaced: list[Path] = []
    try:
        for path, destination in backups.items():
            with destination.open("xb") as handle:
                handle.write(originals[path])
                handle.flush()
                os.fsync(handle.fileno())
            shutil.copystat(path, destination)
            created_backups.append(destination)
        for path, temporary in staged.items():
            os.replace(temporary, path)
            replaced.append(path)
        return list(backups.values())
    except Exception:
        for path in replaced:
            recovery = stage_file(path, originals[path])
            os.replace(recovery, path)
        raise
    finally:
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)


def validate_and_render(
    project_root: Path,
    recompute: Callable[[Path], tuple[list[dict[str, Any]], dict[str, Any]]] | None = None,
) -> tuple[Path, Path, str, str]:
    project_root = project_root.resolve()
    run_dir = project_root / "results" / "sse_fss_v1"
    markdown_path = project_root / "paper" / "manuscript_v1.md"
    tex_path = project_root / "paper" / "manuscript_v1.tex"
    figure_base = project_root / "figures" / "figure3_sse"

    summary = load_json(run_dir / "fss_summary.json")
    fit_summary = load_json(run_dir / "fit_results_summary.json")
    rows = load_fss_csv(run_dir / "fss_observables.csv")
    validate_fss_summary(summary)
    validate_fit_summary(fit_summary, rows, summary)
    require_figures(figure_base)
    raw_rows, raw_summary = (recompute or recompute_fss)(run_dir)
    compare_recomputed(rows, summary, raw_rows, raw_summary)

    if not markdown_path.is_file() or not tex_path.is_file():
        raise BackfillError("缺少 Markdown 或 LaTeX 论文源文件")
    markdown = markdown_path.read_text(encoding="utf-8")
    tex = tex_path.read_text(encoding="utf-8")
    validate_placeholders(markdown, tex)
    md_replacements, tex_replacements = build_replacements(fit_summary)
    return (
        markdown_path,
        tex_path,
        replace_once(markdown, md_replacements),
        replace_once(tex, tex_replacements),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-only", action="store_true", help="只检查，不写入")
    mode.add_argument("--apply", action="store_true", help="备份后正式回填")
    args = parser.parse_args(argv)
    try:
        markdown_path, tex_path, markdown, tex = validate_and_render(args.project_root)
        if not args.apply:
            print("正式 SSE 数据和四组论文占位符均通过检查；未写入文件。")
            return 0
        backups = apply_transaction({markdown_path: markdown, tex_path: tex})
        print("正式 SSE 论文回填完成。备份文件：")
        for path in backups:
            print(f"  {path}")
        return 0
    except (RuntimeError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f"正式 SSE 论文回填被拒绝：{exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
