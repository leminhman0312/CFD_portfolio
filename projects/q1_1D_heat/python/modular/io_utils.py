import os
from typing import Tuple
import numpy as np



# ============================================================================
# Helpers for directories
# ============================================================================

def create_dir() -> None:
    os.makedirs("data", exist_ok=True)
    os.makedirs("plot", exist_ok=True)


def make_dt_tag(dt: float) -> str:
    code = int(round(dt * 100.0))
    return f"{code:03d}"


def make_t_tag(t: float) -> str:
    code = int(round(t * 100.0))
    return f"{code:03d}"


def build_grid(imax: int, dx: float) -> np.ndarray:
    x = np.zeros(imax + 1, dtype=float)
    x[1] = 0.0
    for i in range(1, imax):
        x[i + 1] = x[i] + dx
    return x


def load_block_from_file(path: str, idx: int):
    if not os.path.isfile(path):
        return None, None, False

    block = -1
    x = []
    y = []
    current_t = 0.0

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("# t ="):
                block += 1
                if block > idx:
                    break
                try:
                    current_t = float(s.split("=")[1].split()[0])
                except Exception:
                    current_t = 0.0
                x = []
                y = []
                continue

            if s.startswith("#"):
                continue

            if block == idx:
                parts = s.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))

    if block < idx:
        return None, None, False

    return current_t, (np.array(x), np.array(y)), True


def write_block_T(f, t: float, x: np.ndarray, T: np.ndarray, imax: int) -> None:
    f.write(f"# t = {t:.2f} hr\n")
    f.write("# x\tT\n")
    for i in range(1, imax + 1):
        f.write(f"{x[i]:.6f}\t{T[i]:.6f}\n")
    f.write("\n\n")


def write_block_error(f, t: float, x: np.ndarray, Terr: np.ndarray, imax: int) -> None:
    f.write(f"# t = {t:.2f} hr\n")
    f.write("# x\terr\n")
    for i in range(1, imax + 1):
        f.write(f"{x[i]:.6f}\t{Terr[i]:.6f}\n")
    f.write("\n\n")
