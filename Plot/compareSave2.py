#!/usr/bin/env python3
"""
Compare legacy DATA/DATA_2D/name_dim.mco
with modular Output/Output_2D/itXXXX_name_dim.mco.

Handles different shapes by using normalized coordinates.
The legacy array is interpolated onto the modular grid before subtraction.

Usage:
    python compare_mco.py phi
    python compare_mco.py phi 2001
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII .mco file.

    The first line is assumed to be:
        n1 n2

    The rest of the file is read directly, whatever the shape is.
    """
    data = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(f"Bad header in {path}: expected 'n1 n2', got {header}")

        for line in f:
            row = line.strip().split()
            if row:
                data.append([float(x) for x in row])

    if not data:
        raise ValueError(f"No data found in {path}")

    ncols = len(data[0])

    for i, row in enumerate(data):
        if len(row) != ncols:
            raise ValueError(
                f"{path}: row {i + 1} has {len(row)} values, expected {ncols}"
            )

    return np.array(data, dtype=float)


def resample_to_shape(arr: np.ndarray, target_shape: tuple[int, int]) -> np.ndarray:
    """
    Resample a 2D array onto a new shape using normalized coordinates.

    This uses only numpy, no scipy.
    """
    ny_old, nx_old = arr.shape
    ny_new, nx_new = target_shape

    x_old = np.linspace(0.0, 1.0, nx_old)
    y_old = np.linspace(0.0, 1.0, ny_old)

    x_new = np.linspace(0.0, 1.0, nx_new)
    y_new = np.linspace(0.0, 1.0, ny_new)

    # Interpolate each row in x.
    tmp = np.empty((ny_old, nx_new))

    for j in range(ny_old):
        tmp[j, :] = np.interp(x_new, x_old, arr[j, :])

    # Interpolate each column in y.
    out = np.empty((ny_new, nx_new))

    for i in range(nx_new):
        out[:, i] = np.interp(y_new, y_old, tmp[:, i])

    return out


def axis_labels(dim: str) -> tuple[str, str]:
    if dim == "xy":
        return "normalized x", "normalized y"
    if dim == "xz":
        return "normalized x", "normalized z"
    if dim == "yz":
        return "normalized y", "normalized z"

    return "normalized axis 1", "normalized axis 2"


def compare_one_dim(name: str, dim: str, iteration: int, root_dir: Path) -> None:
    legacy_path = root_dir / "DATA" / "DATA_2D" / f"{name}_{dim}.mco"
    modular_path = (
        root_dir / "Output" / "Output_2D" / f"it{iteration}_{name}_{dim}.mco"
    )

    if not legacy_path.exists():
        raise FileNotFoundError(f"Missing legacy file: {legacy_path}")

    if not modular_path.exists():
        raise FileNotFoundError(f"Missing modular file: {modular_path}")

    arr_legacy = read_mco_ascii(legacy_path)
    arr_modular = read_mco_ascii(modular_path)

    print(f"\n{dim}:")
    print(f"  legacy shape  = {arr_legacy.shape}")
    print(f"  modular shape = {arr_modular.shape}")

    # Put legacy on modular grid before subtraction.
    arr_legacy_interp = resample_to_shape(arr_legacy, arr_modular.shape)
    diff = arr_modular - arr_legacy_interp

    print(f"  comparison shape = {diff.shape}")
    print(f"  max abs diff     = {np.max(np.abs(diff)):.6e}")
    print(f"  mean diff        = {np.mean(diff):.6e}")

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(1, 4, figsize=(24, 5))

    extent = (0.0, 1.0, 0.0, 1.0)

    im1 = axes[0].imshow(
        arr_legacy,
        origin="lower",
        aspect="auto",
        extent=extent,
    )
    plt.colorbar(im1, ax=axes[0])
    axes[0].set_title(f"Legacy\n{legacy_path.name}\nshape={arr_legacy.shape}")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    im2 = axes[1].imshow(
        arr_modular,
        origin="lower",
        aspect="auto",
        extent=extent,
    )
    plt.colorbar(im2, ax=axes[1])
    axes[1].set_title(f"Modular\n{modular_path.name}\nshape={arr_modular.shape}")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(ylabel)

    im3 = axes[2].imshow(
        diff,
        origin="lower",
        aspect="auto",
        extent=extent,
    )
    plt.colorbar(im3, ax=axes[2])
    axes[2].set_title("Modular - legacy interpolated")
    axes[2].set_xlabel(xlabel)
    axes[2].set_ylabel(ylabel)

    # 1D profiles: average along vertical direction.
    avg_legacy = np.mean(arr_legacy, axis=0)
    avg_modular = np.mean(arr_modular, axis=0)

    x_legacy = np.linspace(0.0, 1.0, avg_legacy.size)
    x_modular = np.linspace(0.0, 1.0, avg_modular.size)

    avg_legacy_interp = np.interp(x_modular, x_legacy, avg_legacy)

    axes[3].plot(
        x_modular,
        avg_legacy_interp,
        label="legacy interpolated",
        linewidth=3,
    )
    axes[3].plot(
        x_modular,
        avg_modular,
        label="modular",
        linewidth=3,
    )

    axes[3].set_title(f"Average profiles, {dim}")
    axes[3].set_xlabel(xlabel)
    axes[3].set_ylabel(name)
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    fig.suptitle(f"{name}_{dim} — normalized comparison", fontsize=14)

    plt.tight_layout()

    output_file = root_dir / f"compare_{name}_{dim}.png"
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"  saved: {output_file}")


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage:")
        print("    python compare_mco.py name")
        print("    python compare_mco.py name iteration")
        print("")
        print("Example:")
        print("    python compare_mco.py phi 2001")
        return 1

    name = sys.argv[1]

    if len(sys.argv) >= 3:
        iteration = int(sys.argv[2])
    else:
        iteration = 1001

    root_dir = Path.cwd()

    for dim in ("xy", "xz", "yz"):
        try:
            compare_one_dim(name, dim, iteration, root_dir)
        except Exception as e:
            print(f"Error for dim={dim}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())