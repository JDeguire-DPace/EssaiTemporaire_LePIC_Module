import numpy as np
from pathlib import Path


def read_mco(path: str):
    path = Path(path)
    with path.open("r") as f:
        header = f.readline().split()
        nx, ny = map(int, header[:2])
        nx+=1
        ny+=1
        data = np.loadtxt(f)

    if data.shape != (ny, nx):
        raise ValueError(
            f"{path}: expected shape {(ny, nx)} from header, got {data.shape}"
        )

    return data


def compare_mco(file_a: str, file_b: str, name: str = ""):
    a = read_mco(file_a)
    b = read_mco(file_b)

    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")

    d = b - a

    print(f"\n=== {name or Path(file_a).name} ===")
    print(f"shape         = {a.shape}")
    print(f"max abs diff  = {np.max(np.abs(d)):.16e}")
    print(f"mean abs diff = {np.mean(np.abs(d)):.16e}")
    print(f"rms diff      = {np.sqrt(np.mean(d**2)):.16e}")
    print(f"sum diff      = {np.sum(d):.16e}")
    print(f"sum abs diff  = {np.sum(np.abs(d)):.16e}")


if __name__ == "__main__":
    # Example for phi planes
    compare_mco("DATA/DATA_2D/phi_xy.mco", "Output/Output_2D/phi1_xy.mco", "phi_xy")
    compare_mco("DATA/DATA_2D/phi_xz.mco", "Output/Output_2D/phi1_xz.mco", "phi_xy")
    compare_mco("DATA/DATA_2D/phi_yz.mco", "Output/Output_2D/phi1_yz.mco", "phi_xy")