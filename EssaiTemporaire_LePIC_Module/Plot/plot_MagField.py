#!/usr/bin/env python3
"""
Plot an ASCII .mco plane file written by mod_io_legacy_mco.f90.

Expected format:
  line 1: nx ny
  then ny lines, each with nx floating-point values (space-separated)

Works for xy/xz/yz planes (just a 2D array).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Tkinter file dialog (built-in)
import tkinter as tk
from tkinter import filedialog, messagebox


def pick_file(initial_dir: Path) -> Path | None:
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    file_path = filedialog.askopenfilename(
        title="Select an .mco file to plot",
        initialdir=str(initial_dir),
        filetypes=[
            ("MCO files", "*.mco"),
            ("All files", "*.*"),
        ],
    )
    root.destroy()

    if not file_path:
        return None
    return Path(file_path)


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII mco: first line 'n1 n2', then n2 rows with n1 values.
    Returns array with shape (n2, n1) for plotting (row-major).
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path.name}: expected 2 ints 'n1 n2', got: {header}"
            )

        n1 = int(header[0]) + 1
        n2 = int(header[1]) + 1

        data = []
        for _ in range(n2):
            line = f.readline()
            if not line:
                raise ValueError(
                    f"Unexpected EOF in {path.name}: expected {n2} data rows."
                )
            row = line.strip().split()
            if len(row) != n1:
                raise ValueError(
                    f"Row has {len(row)} values but expected {n1} in {path.name}."
                )
            data.append([float(x) for x in row])

    return np.array(data, dtype=float)


    

def main() -> int:
    # def essai_fit(x,y,z):
    #     x0 = 0.5*128*3.125e-3
    #     y0 = 0.5*96*3.3333e-3
    #     z0 = 0.5*192*3.0208e-3
    #     sig2 = (0.05)**2
    #     eps0 = 8.854187817e-12
    #     A = (1.0/np.sqrt(2.*np.pi*sig2))*(1.6022e-19 * 1e12/eps0)
    #     return A * np.exp(-((x-x0)**2 + (y-y0)**2 + (z-z0)**2)/(2*sig2))
    which_directory = input("DATA or Output? \n")
    Dir2D = which_directory+"_2D"
    script_dir = Path(__file__).resolve().parent
    default_dir = (script_dir.parent / which_directory / Dir2D)
    initial_dir = default_dir if default_dir.is_dir() else script_dir.parent

    # xo = np.linspace(0, 0.4, 128)
    # yo = np.linspace(0, 0.32, 96)
    # zo = np.linspace(0, 0.56, 192)

    # phi = essai_fit(*np.meshgrid(xo, yo, zo, indexing="ij"))

    # plt.figure()
    # plt.plot(xo, phi[:, 48, 96])
    #im = plt.imshow(phi[:, 48, :], origin="lower", aspect="auto")
    #plt.colorbar(im, label="Phi")
    #plt.title("Test Gaussian Source (Log scale)")
    #plt.xlabel("Index 1")
    #plt.ylabel("Index 2")
    #plt.tight_layout()
    

    # Discover files starting with 'B' in the DATA_2D directory.
    files = []
    if default_dir.is_dir():
        files = sorted(default_dir.glob("B*"))

    # If no B* files found, fall back to asking the user for a single file.
    if not files:
        path = pick_file(initial_dir)
        if path is None:
            print("No file selected. Exiting.")
            return 0
        files = [path]
    else:
        files = files[:9]

    arrays: list[np.ndarray] = []
    names: list[str] = []
    for p in files:
        try:
            arrays.append(read_mco_ascii(p))
            names.append(p.name)
        except Exception as e:
            print(f"Failed to read {p}: {e}", file=sys.stderr)

    if not arrays:
        print("No valid data to plot. Exiting.")
        return 1

    # Use a consistent color scale across all subplots
    vmin = min(a.min() for a in arrays)
    vmax = max(a.max() for a in arrays)

    n = len(arrays)
    cols = 3
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    if isinstance(axes, np.ndarray):
        axes_flat = axes.flatten()
    else:
        axes_flat = [axes]

    im = None
    for i, ax in enumerate(axes_flat):
        if i < n:
            im = ax.imshow(arrays[i], origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
            ax.set_title(names[i])
            ax.set_xlabel("x cell")
            ax.set_ylabel("z cell")
        else:
            ax.axis("off")

    # Shared colorbar for all subplots
    if im is not None:
        fig.colorbar(im, ax=axes_flat[:n], orientation="vertical", fraction=0.02, pad=0.04, label="B")

    plt.tight_layout()
    plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
