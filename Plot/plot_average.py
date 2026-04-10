#!/usr/bin/env python3
"""
Average multiple ASCII .mco files.

Expected files:
  it11_phi_xz.mco, it21_phi_xz.mco, ..., it191_phi_xz.mco

Output:
  - Displays averaged field
  - Saves avg_phi_xz.mco
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------
# READ MCO
# ---------------------------------
def read_mco_ascii(path: Path) -> np.ndarray:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(f"Bad header in {path.name}")

        n1 = int(header[0]) + 1
        n2 = int(header[1]) + 1

        data = []
        for _ in range(n2):
            row = f.readline().split()
            if len(row) != n1:
                raise ValueError(f"Bad row length in {path.name}")
            data.append([float(x) for x in row])

    return np.array(data)


# ---------------------------------
# WRITE MCO
# ---------------------------------
def write_mco_ascii(path: Path, arr: np.ndarray):
    n2, n1 = arr.shape

    with path.open("w") as f:
        f.write(f"{n1-1} {n2-1}\n")
        for j in range(n2):
            f.write(" ".join(f"{v:.6e}" for v in arr[j, :]) + "\n")


# ---------------------------------
# MAIN
# ---------------------------------
def main():

    # 🔧 CHANGE THIS IF NEEDED
    folder = Path(__file__).resolve().parent

    print(f"Looking in: {folder}\n")

    # Collect files
    files = []
    for i in range(1, 30):
        f = folder / f"../Output/Output_2D/it{i}1_n1_xz.mco"
        if f.exists():
            files.append(f)
        else:
            print(f"Missing: {f.name}")

    if not files:
        raise RuntimeError("No files found.")

    print("\nUsing files:")
    for f in files:
        print(f"  {f.name}")

    # Read arrays
    arrays = []
    for f in files:
        try:
            arr = read_mco_ascii(f)
            arrays.append(arr)
            print(f"Read OK: {f.name}")
        except Exception as e:
            print(f"Skipping {f.name}: {e}")

    if not arrays:
        raise RuntimeError("No valid files could be read.")

    # Check shapes
    shape0 = arrays[0].shape
    for arr in arrays:
        if arr.shape != shape0:
            raise ValueError("Shape mismatch between files.")

    # ---------------------------------
    # AVERAGE
    # ---------------------------------
    avg = np.mean(arrays, axis=0)

    print(f"\nAveraged {len(arrays)} files.")

    # ---------------------------------
    # SAVE RESULT
    # ---------------------------------
    out_file = folder / "avg_phi_xz.mco"
    write_mco_ascii(out_file, avg)
    print(f"Saved: {out_file}")

    # ---------------------------------
    # PLOT
    # ---------------------------------
    plt.figure(figsize=(7, 5))

    im = plt.imshow(avg, origin="lower", aspect="auto")
    plt.colorbar(im, label=r"$\phi$ (V)")
    plt.title("Averaged Potential")
    plt.xlabel("x cell")
    plt.ylabel("z cell")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()