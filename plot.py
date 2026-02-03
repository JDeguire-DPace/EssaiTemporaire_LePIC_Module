#!/usr/bin/env python3
"""
read_bcnd_mco_and_plot.py

- Reads a bcnd-style ASCII file:
    nx ny
    <ny+1 lines, each with nx+1 integers>

- Loads it into a NumPy array of shape (ny+1, nx+1)
- Optionally plots it with matplotlib, allowing cmap and color scale control.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def read_bcnd_file(path):
    """Read bcnd-style text file -> 2D numpy array."""
    with open(path, "r") as f:
        # --- Header: nx ny ---
        header = f.readline().strip().split()
        if len(header) < 2:
            raise ValueError(f"Header line must contain at least 2 integers, got: {header}")

        nx = int(header[0])
        ny = int(header[1])

        ncols = nx + 1
        nrows = ny + 1
        total_values_expected = nrows * ncols

        # --- Read all remaining ints ---
        values = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            values.extend(line.split())

        if len(values) != total_values_expected:
            raise ValueError(
                f"Expected {total_values_expected} values after header, "
                f"but found {len(values)}"
            )

        data = np.array(values, dtype=float).reshape((nrows, ncols))
        return data


def main():
    parser = argparse.ArgumentParser(description="Read and plot a bcnd-style ASCII array file.")
    parser.add_argument("input", help="Path to the text file (e.g. bcnd_xy.mco)")
    parser.add_argument("--save-npy", metavar="OUT.npy",
                        help="Optional: save the array as a NumPy .npy file")
    parser.add_argument("--show-plot", action="store_true",
                        help="Show a matplotlib plot of the array")

    # Colormap and scale options
    parser.add_argument("--cmap", default="seismic",
                        help="Matplotlib colormap name (default: seismic)")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Minimum value for color scale (default: automatic)")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Maximum value for color scale (default: automatic)")

    args = parser.parse_args()

    arr = read_bcnd_file(args.input)
    print(f"Loaded array from {args.input}")
    print(f"Shape: {arr.shape} (rows, cols)")
    print("Example values:")
    print(arr[:5, :10])
    plt.figure(figsize=(6, 5))
    im = plt.imshow(arr, origin="lower", cmap='plasma', norm=LogNorm())
    plt.colorbar()
    plt.show()    


if __name__ == "__main__":
    main()
