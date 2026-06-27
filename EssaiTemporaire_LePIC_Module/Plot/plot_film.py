#!/usr/bin/env python3
"""
Animate ASCII .mco xz plane files written by mod_io_legacy_mco.f90.

Expected format:
  line 1: n1 n2
  then n2+1 lines, each with n1+1 floating-point values

Expected filenames:
  it11_phi_xz.mco, it21_phi_xz.mco, ..., it191_phi_xz.mco
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII mco: first line 'n1 n2', then n2+1 rows with n1+1 values.
    Returns array with shape (n2+1, n1+1).
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path.name}: expected 2 ints 'n1 n2', got: {header}"
            )

        n1 = int(header[0]) 
        n2 = int(header[1]) 

        data = []
        for j in range(n2):
            line = f.readline()
            if not line:
                raise ValueError(
                    f"Unexpected EOF in {path.name}: expected {n2} data rows."
                )

            row = line.strip().split()
            if len(row) != n1:
                raise ValueError(
                    f"Row {j+1} in {path.name} has {len(row)} values but expected {n1}."
                )

            data.append([float(x) for x in row])

    return np.array(data, dtype=float)


def main() -> int:
    folder = Path(__file__).resolve().parent

    print(f"Looking for files in:\n  {folder}\n")

    file_list = [folder / f"../Output/Output_2D/it{i}_n1_xz.mco" for i in range(1, 119002, 1000)]

    existing_files = []
    for f in file_list:
        if f.exists():
            existing_files.append(f)
        else:
            print(f"Missing: {f.name}")

    if not existing_files:
        raise FileNotFoundError(
            "No matching files found.\n"
            "Expected names like: it11_phi_xz.mco, it21_phi_xz.mco, ..., it191_phi_xz.mco\n"
            f"Search folder: {folder}"
        )

    print("\nFiles used in animation:")
    for f in existing_files:
        print(f"  {f.name}")

    data_list = []
    used_files = []

    for f in existing_files:
        try:
            arr = read_mco_ascii(f)
            data_list.append(arr)
            used_files.append(f)
            print(f"Read OK: {f.name}   shape={arr.shape}")
        except Exception as e:
            print(f"Skipping {f.name}: {e}")

    if not data_list:
        raise RuntimeError(
            "Files were found, but none could be read successfully as ASCII .mco files."
        )

    shape0 = data_list[0].shape
    for f, arr in zip(used_files, data_list):
        if arr.shape != shape0:
            raise ValueError(
                f"Shape mismatch: {f.name} has shape {arr.shape}, expected {shape0}"
            )

    frame_duration_s = 0.16
    fps = 1.0 / frame_duration_s

    fig, ax = plt.subplots(figsize=(8, 6))

    # Initial limits from first frame
    vmin0 = np.min(data_list[0])
    vmax0 = np.max(data_list[0])

    im = ax.imshow(
        data_list[0],
        origin="lower",
        aspect="auto",
        vmin=vmin0,
        vmax=vmax0,
        cmap="seismic",
    )

    cbar = plt.colorbar(im, ax=ax, label=r"$\phi$ (V)")
    ax.set_xlabel("x cell")
    ax.set_ylabel("z cell")
    title = ax.set_title(used_files[0].name)

    def update(frame: int):
        data = data_list[frame]
        vmin = np.min(data)
        vmax = np.max(data)

        # Avoid identical min/max which breaks the color scale
        if np.isclose(vmin, vmax):
            delta = 1e-12 if np.isclose(vmin, 0.0) else 1e-12 * abs(vmin)
            vmin -= delta
            vmax += delta

        im.set_data(data)
        im.set_clim(vmin, vmax)
        cbar.update_normal(im)

        title.set_text(
            f"{used_files[frame].name}   "
            f"(min={np.min(data):.3e}, max={np.max(data):.3e})"
        )
        return [im, title]

    ani = FuncAnimation(
        fig,
        update,
        frames=len(data_list),
        interval=frame_duration_s * 1000,
        blit=False,
        repeat=True,
    )

    gif_path = folder / "phi_xz_animation.gif"
    print(f"\nSaving GIF to: {gif_path}")
    ani.save(gif_path, writer=PillowWriter(fps=fps))

    plt.tight_layout()
    plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())