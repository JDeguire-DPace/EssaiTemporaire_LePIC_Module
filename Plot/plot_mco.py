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

        n1 = int(header[0])+1
        n2 = int(header[1])+1

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

    arr = np.array(data, dtype=float)  # shape (n2, n1)
    return arr


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_dir = (script_dir.parent / "DATA" / "DATA_2D")
    initial_dir = default_dir if default_dir.is_dir() else script_dir.parent

    path = pick_file(initial_dir)
    if path is None:
        print("No file selected. Exiting.")
        return 0

    try:
        arr = read_mco_ascii(path)
    except Exception as e:
        # show GUI error and also print to terminal
        try:
            tk.Tk().withdraw()
            messagebox.showerror("Failed to read .mco", str(e))
        except Exception:
            pass
        print(f"Error: {e}", file=sys.stderr)
        return 1

    plt.figure()
    plt.imshow(arr, origin="lower", aspect="auto")
    plt.colorbar(label=path.stem)
    plt.title(path.name)
    plt.xlabel("Index 1")
    plt.ylabel("Index 2")
    plt.tight_layout()
    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
