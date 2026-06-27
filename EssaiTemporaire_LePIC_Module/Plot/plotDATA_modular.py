#!/usr/bin/env python3
"""
Plot modular Output/Output_2D/itXXXX_name_dim.mco files.

Usage:
    python plot_mco_modular.py phi
    python plot_mco_modular.py phi 2001

This plots:
    Output/Output_2D/itXXXX_name_xy.mco
    Output/Output_2D/itXXXX_name_xz.mco
    Output/Output_2D/itXXXX_name_yz.mco

Figures are saved in the current directory.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII .mco file.

    First line:
        n1 n2

    Remaining lines:
        data values
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


def axis_labels(dim: str) -> tuple[str, str]:
    if dim == "xy":
        return "x index", "y index"
    if dim == "xz":
        return "x index", "z index"
    if dim == "yz":
        return "y index", "z index"

    return "axis 1 index", "axis 2 index"


def plot_one_dim(name: str, dim: str, iteration: int, root_dir: Path) -> None:
    modular_path = (
        root_dir / "Output" / "Output_2D" / f"it{iteration}_{name}_{dim}.mco"
    )

    if not modular_path.exists():
        raise FileNotFoundError(f"Missing modular file: {modular_path}")

    arr = read_mco_ascii(modular_path)

    print(f"\n{dim}:")
    print(f"  file  = {modular_path}")
    print(f"  shape = {arr.shape}")
    print(f"  min   = {np.min(arr):.6e}")
    print(f"  max   = {np.max(arr):.6e}")
    print(f"  mean  = {np.mean(arr):.6e}")

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))


    im = axes[0].imshow(
        arr,
        origin="lower",
        aspect="auto",
    )

    plt.colorbar(im, ax=axes[0])
    axes[0].set_title(f"{modular_path.name}\nshape={arr.shape}")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    # 1D profile: average along vertical direction.
    avg_profile = np.mean(arr, axis=0)
    x = np.arange(avg_profile.size)

    axes[1].plot(
        x,
        avg_profile,
        linewidth=3,
    )

    axes[1].set_title(f"Average profile, {dim}")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(name)
    axes[1].grid(True, alpha=0.3)

    fig.suptitle(f"{name}_{dim} modular output", fontsize=14)

    plt.tight_layout()

    output_file = root_dir / f"plot_{name}_{dim}.png"
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"  saved: {output_file}")


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage:")
        print("    python plot_mco_modular.py name")
        print("    python plot_mco_modular.py name iteration")
        print("")
        print("Example:")
        print("    python plot_mco_modular.py phi 2001")
        return 1

    name = sys.argv[1]

    if len(sys.argv) >= 3:
        iteration = int(sys.argv[2])
    else:
        iteration = 118001

    root_dir = Path.cwd()

    for dim in ("xy", "xz", "yz"):
        try:
            plot_one_dim(name, dim, iteration, root_dir)
        except Exception as e:
            print(f"Error for dim={dim}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())