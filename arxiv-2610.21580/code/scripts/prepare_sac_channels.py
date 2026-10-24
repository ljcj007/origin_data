#!/usr/bin/env python3
"""生成 super-honeycomb 动力学 SSE 的动量与探针通道。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_graph(path: Path) -> tuple[int, list[tuple[int, int, int, int]], np.ndarray]:
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    nbasis, nedge, dimension = map(int, lines[0])
    if dimension != 3 or len(lines) != nedge + 3:
        raise ValueError("格点图格式错误")
    edges = []
    for row in lines[1 : nedge + 1]:
        site1, site2, dx, dy, _dz, _kind = map(int, row)
        edges.append((site1 - 1, site2 - 1, dx, dy))
    eta = np.asarray([int(value) for value in lines[nedge + 1]], dtype=float)
    if eta.shape != (nbasis,) or not np.all(np.isin(eta, (-1.0, 1.0))):
        raise ValueError("需要带二分子格标记的格点图")
    return nbasis, edges, eta


def read_positions(path: Path, nbasis: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    atoms = sorted(payload["atoms"], key=lambda item: item["site"])
    positions = np.asarray([item["fractional"][:2] for item in atoms], dtype=float)
    if positions.shape != (nbasis, 2):
        raise ValueError("原胞位置数与格点图不一致")
    return positions


def reciprocal_path(size: int) -> list[tuple[int, int, str]]:
    if size < 6 or size % 6:
        raise ValueError("Γ-M-K-Γ 离散路径要求 L 是 6 的倍数")
    points: list[tuple[int, int, str]] = []
    for mx in range(size // 2 + 1):
        points.append((mx, 0, "Γ" if mx == 0 else ("M" if mx == size // 2 else "path")))
    for step in range(1, size // 6 + 1):
        mx = size // 2 - step
        my = 2 * step
        points.append((mx, my, "K" if step == size // 6 else "path"))
    for step in range(1, size // 3 + 1):
        mx = size // 3 - step
        my = size // 3 - step
        points.append((mx, my, "Γ" if step == size // 3 else "path"))
    return points


def coupling_matrix(
    nbasis: int, edges: list[tuple[int, int, int, int]], eta: np.ndarray, kx: float, ky: float
) -> tuple[np.ndarray, list[int], list[int]]:
    plus = [index for index in range(nbasis) if eta[index] > 0]
    minus = [index for index in range(nbasis) if eta[index] < 0]
    plus_index = {site: index for index, site in enumerate(plus)}
    minus_index = {site: index for index, site in enumerate(minus)}
    matrix = np.zeros((len(plus), len(minus)), dtype=np.complex128)
    for site1, site2, dx, dy in edges:
        if eta[site1] > 0:
            matrix[plus_index[site1], minus_index[site2]] += np.exp(1j * (kx * dx + ky * dy))
        else:
            matrix[plus_index[site2], minus_index[site1]] += np.exp(-1j * (kx * dx + ky * dy))
    return matrix, plus, minus


def phase_vector(kx: float, ky: float, positions: np.ndarray) -> np.ndarray:
    return np.exp(1j * (kx * positions[:, 0] + ky * positions[:, 1]))


def channels_at_q(
    nbasis: int,
    edges: list[tuple[int, int, int, int]],
    eta: np.ndarray,
    positions: np.ndarray,
    kx: float,
    ky: float,
    include_sites: bool,
) -> list[tuple[str, np.ndarray]]:
    phase = phase_vector(kx, ky, positions)
    hub = np.where(eta > 0, phase, 0.0)
    edge = np.where(eta < 0, phase, 0.0)
    matrix, _plus, minus = coupling_matrix(nbasis, edges, eta, kx, ky)
    _u, _s, vh = np.linalg.svd(matrix, full_matrices=True)
    null_vector = vh[-1].conjugate()
    pivot = int(np.argmax(np.abs(null_vector)))
    null_vector *= np.exp(-1j * np.angle(null_vector[pivot]))
    cls = np.zeros(nbasis, dtype=np.complex128)
    cls[minus] = null_vector
    channels = [
        ("total", phase),
        ("staggered", eta * phase),
        ("hub", hub),
        ("edge", edge),
        ("cls", cls),
    ]
    if include_sites:
        for site in range(nbasis):
            vector = np.zeros(nbasis, dtype=np.complex128)
            vector[site] = phase[site]
            channels.append((f"site{site + 1}", vector))
    return channels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=Path("graphs/super_honeycomb_poscar.dat"))
    parser.add_argument(
        "--metadata", type=Path, default=Path("results/poscar2/lattice_graph_nn.dat.json")
    )
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--without-sites", action="store_true")
    args = parser.parse_args()

    graph = args.graph.resolve()
    metadata = args.metadata.resolve()
    nbasis, edges, eta = read_graph(graph)
    positions = read_positions(metadata, nbasis)
    path = reciprocal_path(args.size)
    rows = []
    labels: list[str] | None = None
    for q_index, (mx, my, point_label) in enumerate(path, start=1):
        kx = 2.0 * math.pi * mx / args.size
        ky = 2.0 * math.pi * my / args.size
        channels = channels_at_q(
            nbasis, edges, eta, positions, kx, ky, not args.without_sites
        )
        current_labels = [name for name, _vector in channels]
        if labels is None:
            labels = current_labels
        elif current_labels != labels:
            raise RuntimeError("通道顺序发生变化")
        rows.append((q_index, mx, my, point_label, kx, ky, channels))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(rows)} {len(labels or [])} {nbasis}\n")
        for q_index, _mx, _my, _point_label, kx, ky, channels in rows:
            handle.write(f"{q_index} {kx:.17g} {ky:.17g}\n")
            for name, vector in channels:
                values = " ".join(f"{value.real:.17g} {value.imag:.17g}" for value in vector)
                handle.write(f"{name} {values}\n")

    manifest = args.manifest or args.output.with_suffix(".json")
    payload = {
        "format": "sse-dynamics-channels-v1",
        "size": args.size,
        "basis_count": nbasis,
        "path": "Gamma-M-K-Gamma",
        "q_points": [
            {"q_index": q, "mx": mx, "my": my, "label": label, "kx": kx, "ky": ky}
            for q, mx, my, label, kx, ky, _channels in rows
        ],
        "channels": labels,
        "graph": str(graph),
        "graph_sha256": sha256(graph),
        "metadata": str(metadata),
        "metadata_sha256": sha256(metadata),
        "channel_file": str(args.output.resolve()),
        "channel_file_sha256": sha256(args.output),
    }
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"动量点：{len(rows)}，通道：{len(labels or [])}")
    print(f"通道文件：{args.output}")
    print(f"清单：{manifest}")


if __name__ == "__main__":
    main()
