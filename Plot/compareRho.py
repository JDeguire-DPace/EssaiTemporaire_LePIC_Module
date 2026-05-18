#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Settings
# ============================================================

ROOT = Path(".").resolve()
LEGACY_DIR = ROOT / "DATA" / "DATA_2D"
MODULAR_DIR = ROOT / "Output" / "Output_2D"

ITERATION = 1001
DIMS = ("xy", "xz", "yz")

E_CHARGE = 1.602176634e-19

# ============================================================
# Reader
# ============================================================

def read_mco_ascii(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    with path.open("r") as f:
        n1, n2 = map(int, f.readline().split())
        n1 += 1
        n2 += 1

        data = []
        for _ in range(n2):
            row = list(map(float, f.readline().split()))
            data.append(row)

    return np.array(data)


# ============================================================
# Physics
# ============================================================

def build_rho(n1_path, n2_path):
    n1 = read_mco_ascii(n1_path)
    n2 = read_mco_ascii(n2_path)

    return E_CHARGE * (n2 - n1)


# ============================================================
# Plot
# ============================================================

def compare_dim(dim):
    legacy_n1 = LEGACY_DIR / f"n1_{dim}.mco"
    legacy_n2 = LEGACY_DIR / f"n2_{dim}.mco"

    modular_n1 = MODULAR_DIR / f"it{ITERATION}_n1_{dim}.mco"
    modular_n2 = MODULAR_DIR / f"it{ITERATION}_n2_{dim}.mco"

    print(f"\n=== {dim} ===")

    rho_L = build_rho(legacy_n1, legacy_n2)
    rho_M = build_rho(modular_n1, modular_n2)

    diff = rho_M - rho_L

    avg_L = np.mean(rho_L, axis=0)
    avg_M = np.mean(rho_M, axis=0)

    x = np.arange(rho_L.shape[1])

    vmax = np.max(np.abs([rho_L, rho_M]))
    dmax = np.max(np.abs(diff))

    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    im0 = ax[0,0].imshow(rho_L, origin="lower", vmin=-vmax, vmax=vmax)
    ax[0,0].set_title("Legacy")
    plt.colorbar(im0, ax=ax[0,0])

    im1 = ax[0,1].imshow(rho_M, origin="lower", vmin=-vmax, vmax=vmax)
    ax[0,1].set_title("Modular")
    plt.colorbar(im1, ax=ax[0,1])

    im2 = ax[1,0].imshow(diff, origin="lower", vmin=-dmax, vmax=dmax)
    ax[1,0].set_title("Difference")
    plt.colorbar(im2, ax=ax[1,0])

    ax[1,1].plot(x, avg_L, label="Legacy")
    ax[1,1].plot(x, avg_M, label="Modular")
    ax[1,1].set_title("Average profile")
    ax[1,1].legend()
    ax[1,1].grid()

    plt.tight_layout()

    out = f"compare_rho_{dim}.png"
    plt.savefig(out)
    print(f"Saved: {out}")

    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    print("Running rho comparison...")

    for dim in DIMS:
        compare_dim(dim)

    print("Done.")


if __name__ == "__main__":
    main()