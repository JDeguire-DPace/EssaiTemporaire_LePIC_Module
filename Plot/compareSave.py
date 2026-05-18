#!/usr/bin/env python3
"""
Compare legacy DATA/DATA_2D/name_dim.mco
with modular Output/Output_2D/it2001_name_dim.mco.

Usage:
    python compare_mco.py phi

This compares:
    DATA/DATA_2D/phi_xy.mco
    DATA/DATA_2D/phi_xz.mco
    DATA/DATA_2D/phi_yz.mco

with:
    Output/Output_2D/it2001_phi_xy.mco
    Output/Output_2D/it2001_phi_xz.mco
    Output/Output_2D/it2001_phi_yz.mco

Figures are saved in the current directory.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII mco:
        first line: n1 n2
        then n2+1 rows with n1+1 values

    Returns array with shape (n2+1, n1+1).
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path}: expected 'n1 n2', got {header}"
            )

        n1 = int(header[0]) + 1
        n2 = int(header[1]) + 1

        data = []

        for _ in range(n2):
            line = f.readline()

            if not line:
                raise ValueError(
                    f"Unexpected EOF in {path}: expected {n2} rows"
                )

            row = line.strip().split()

            if len(row) != n1:
                raise ValueError(
                    f"{path}: row has {len(row)} values, expected {n1}"
                )

            data.append([float(x) for x in row])

    return np.array(data, dtype=float)


def axis_labels(dim: str) -> tuple[str, str]:
    if dim == "xy":
        return "x cell", "y cell"
    if dim == "xz":
        return "x cell", "z cell"
    if dim == "yz":
        return "y cell", "z cell"

    return "cell 1", "cell 2"


def compare_one_dim(name: str, dim: str, iteration: int, root_dir: Path) -> None:
    legacy_path = root_dir / "DATA" / "DATA_2D" / f"{name}_{dim}.mco"
    modular_path = root_dir / "Output" / "Output_2D" / f"it{iteration}_{name}_{dim}.mco"

    if not legacy_path.exists():
        raise FileNotFoundError(f"Missing legacy file: {legacy_path}")

    if not modular_path.exists():
        raise FileNotFoundError(f"Missing modular file: {modular_path}")

    arr1 = read_mco_ascii(legacy_path)
    arr2 = read_mco_ascii(modular_path)

    if arr1.shape != arr2.shape:
        print(
            f"Warning: shape mismatch for {dim}: "
            f"legacy {arr1.shape}, modular {arr2.shape}"
        )

    avg1 = np.mean(arr1, axis=0)
    avg2 = np.mean(arr2, axis=0)

    x1 = np.arange(arr1.shape[1])
    x2 = np.arange(arr2.shape[1])

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    im1 = axes[0].imshow(arr1, origin="lower", aspect="auto")
    plt.colorbar(im1, ax=axes[0])
    axes[0].set_title(legacy_path.name)
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    im2 = axes[1].imshow(arr2, origin="lower", aspect="auto")
    plt.colorbar(im2, ax=axes[1])
    axes[1].set_title(modular_path.name)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(ylabel)

    axes[2].plot(x1, avg1, label="legacy", linewidth=3)
    axes[2].plot(x2, avg2, label="modular", linewidth=3)
    axes[2].set_title(f"Average profiles, {dim}")
    axes[2].set_xlabel(xlabel)
    axes[2].set_ylabel(name)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(f"{name}_{dim}", fontsize=14)

    plt.tight_layout()

    output_file = root_dir / f"compare_{name}_{dim}.png"
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"Saved: {output_file}")


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage:")
        print("    python compare_mco.py name")
        print("")
        print("Example:")
        print("    python compare_mco.py phi")
        return 1

    name = sys.argv[1]

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
