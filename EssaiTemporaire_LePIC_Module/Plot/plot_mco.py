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
    

    #path1 = Path(rf"/mnt/c/Users/Jasmin Deguire/Desktop/LePIC_3D/DATA/DATA_2D/n1_xz.mco")
    path1 = Path(rf"/mnt/c/Users/Jasmin Deguire/Desktop/LePIC_3D/Output/Output_2D/it2001_n1_xz.mco")
    arr1 = read_mco_ascii(path1)
    #path2 = Path(rf"/mnt/c/Users/Jasmin Deguire/Desktop/LePIC_3D/DATA/DATA_2D/n2_xz.mco")
    path2 = Path(rf"/mnt/c/Users/Jasmin Deguire/Desktop/LePIC_3D/Output/Output_2D/it2001_n2_xz.mco")
    arr2 = read_mco_ascii(path2)

    arr4 = 1.602e-19*(arr2-arr1)

    path3 = Path(rf"/mnt/c/Users/Jasmin Deguire/Desktop/LePIC_3D/Output/Output_2D/it2001_phi_xz.mco")
    arr3 = read_mco_ascii(path3)
    # Create a single figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- 2D plot (left) ---
    im = ax1.imshow(
        arr4,
        origin="lower",
        aspect="auto",
        #norm=LogNorm(vmin=vmin, vmax=vmax),
    )
    plt.colorbar(im, ax=ax1, label=rf"$\rho$ (C/m$^3$)")
    ax1.set_xlabel("x cell")
    ax1.set_ylabel("z cell")

    # --- 1D slice plot (right) ---
    average = np.zeros(arr1.shape[1])
    for i in range(0, arr1.shape[0]):
        average += arr4[i,:]/arr1.shape[0]
    ax2.plot(np.arange(arr1.shape[1]), average, label="LePIC potential", color="skyblue", linewidth=4)
    ax2.set_xlabel("x cell")
    ax2.set_ylabel(rf"$\phi$ (V)")

    x = np.linspace(0,0.1,128)
    eps0 = 8.854e-12

    # -----------------------------
    # Quadratique
    # -----------------------------
    # rho_0 = 1e-7
    # L = 0.1
    # phi_x = 20 - 1e-7*x*(L - x)/(2*eps0)
    
    # -----------------------------
    # Gaussien
    # -----------------------------
    # amplitude = 1e-12
    # x0 = 0.05
    # sig = 0.05/4

    # phi_x = (amplitude/np.sqrt(2*np.pi*sig**2))*np.exp(-((x)-x0)**2/(2*sig**2))/(eps0)
    # ax2.plot(x*(128/0.1),phi_x)

    # -----------------------------
    # Sinus
    # -----------------------------
    # amplitude = 1e-6
    L = 0.10

    phi_x = np.sin(4*np.pi*(x)/L)
    #ax2.plot(x*(128/0.1),phi_x)
    # ax2.plot(x*(128/0.1), phi_x, label="Analytical", color="red", linestyle=(0, (5, 10)), linewidth=4)
    # ax2.legend()
    plt.tight_layout()
    plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
