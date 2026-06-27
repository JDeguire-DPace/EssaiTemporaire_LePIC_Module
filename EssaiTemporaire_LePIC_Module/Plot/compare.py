#!/usr/bin/env python3
"""
Compare two ASCII .mco plane files.

Features:
- Select two files (even from different directories)
- OR pass them via CLI: python plot.py file1.mco file2.mco
- Plot:
    • 2D map of file 1
    • 2D map of file 2
    • 1D averaged profiles on same plot
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox


# -----------------------------
# File selection
# -----------------------------
def pick_two_files(initial_dir: Path) -> tuple[Path, Path] | None:
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    # First file
    file1 = filedialog.askopenfilename(
        title="Select FIRST .mco file",
        initialdir=str(initial_dir),
        filetypes=[("MCO files", "*.mco"), ("All files", "*.*")]
    )
    if not file1:
        root.destroy()
        return None

    # Second file (can be elsewhere)
    file2 = filedialog.askopenfilename(
        title="Select SECOND .mco file",
        initialdir=str(Path(file1).parent),
        filetypes=[("MCO files", "*.mco"), ("All files", "*.*")]
    )
    root.destroy()

    if not file2:
        return None

    return Path(file1), Path(file2)


# -----------------------------
# MCO reader
# -----------------------------
def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII mco: first line 'n1 n2',
    then n2+1 rows with n1+1 values.

    Returns array shape (n2+1, n1+1)
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path.name}: expected 'n1 n2', got {header}"
            )

        n1 = int(header[0]) + 1
        n2 = int(header[1]) + 1

        data = []
        for _ in range(n2):
            line = f.readline()
            if not line:
                raise ValueError(
                    f"Unexpected EOF in {path.name}: expected {n2} rows"
                )

            row = line.strip().split()
            if len(row) != n1:
                raise ValueError(
                    f"{path.name}: row has {len(row)} values, expected {n1}"
                )

            data.append([float(x) for x in row])

    return np.array(data, dtype=float)


# -----------------------------
# Main
# -----------------------------
def main() -> int:
    script_dir = Path(__file__).resolve().parent
    initial_dir = script_dir.parent if script_dir.parent.is_dir() else script_dir

    # --- CLI mode ---
    if len(sys.argv) == 3:
        path1 = Path(sys.argv[1])
        path2 = Path(sys.argv[2])
    else:
        result = pick_two_files(initial_dir)
        if result is None:
            print("File selection cancelled.")
            return 0
        path1, path2 = result

    # --- Read files ---
    try:
        arr1 = read_mco_ascii(path1)
        arr2 = read_mco_ascii(path2)
    except Exception as e:
        try:
            tk.Tk().withdraw()
            messagebox.showerror("Failed to read .mco", str(e))
        except Exception:
            pass
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # --- 1D averages (z-averaged) ---
    avg1 = np.mean(arr1, axis=0)
    avg2 = np.mean(arr2, axis=0)

    x1 = np.arange(arr1.shape[1])
    x2 = np.arange(arr2.shape[1])

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # --- 2D plot 1 ---
    im1 = axes[0].imshow(arr1, origin="lower", aspect="auto")
    plt.colorbar(im1, ax=axes[0], label=r"$\phi$ (V)")
    axes[0].set_title(path1.name)
    axes[0].set_xlabel("x cell")
    axes[0].set_ylabel("z cell")

    # --- 2D plot 2 ---
    im2 = axes[1].imshow(arr2, origin="lower", aspect="auto")
    plt.colorbar(im2, ax=axes[1], label=r"$\phi$ (V)")
    axes[1].set_title(path2.name)
    axes[1].set_xlabel("x cell")
    axes[1].set_ylabel("z cell")

    # --- 1D comparison ---
    axes[2].plot(x1, avg1, label=path1.name, linewidth=3)
    axes[2].plot(x2, avg2, label=path2.name, linewidth=3)
    axes[2].set_title("Average profiles (z-averaged)")
    axes[2].set_xlabel("x cell")
    axes[2].set_ylabel(r"$\phi$ (V)")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
