from io_utils import create_dir
from driver import sim, convergence_study
from post import plot
  

# ============================================================================
# Main, mirrors the C++ main
# ============================================================================

def main():
    imax = 21
    dx = 0.05
    dt = 0.01
    alpha = 0.1
    t0 = 100.0
    tb = 300.0

    t_target = 0.4
    F = alpha * dt / (dx * dx)

    create_dir()
    sim(dx, dt, imax, t0, tb, F)
    convergence_study(dx, imax, alpha, t0, tb, t_target)
    plot(dt, t_target)

    print("Done")


if __name__ == "__main__":
    main()
