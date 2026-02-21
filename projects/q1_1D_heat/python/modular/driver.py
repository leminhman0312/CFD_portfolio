import os
import numpy as np
import matplotlib.pyplot as plt

from io_utils import create_dir, make_dt_tag, build_grid, write_block_T, write_block_error, make_t_tag

from schemes import FTCS_explicit, DufortFrankel, FTCS_implicit,CrankNicolson

from exact import exact_solution_profile

from norms import err_L2, err_Linf



# ============================================================================
# sim: write profile and error files at t = 0.0 0.1 0.2 0.3 0.4
# ============================================================================

def sim(dx: float, dt: float, imax: int, t0: float, tb: float, F: float) -> None:
    tag = make_dt_tag(dt)

    fexp = os.path.join("data", f"ftcs_explicit_{tag}.txt")
    fdf  = os.path.join("data", f"dufort_{tag}.txt")
    fimp = os.path.join("data", f"ftcs_implicit_{tag}.txt")
    fcn  = os.path.join("data", f"cn_{tag}.txt")
    fex  = os.path.join("data", f"exact_{tag}.txt")

    ferr_exp = os.path.join("data", f"error_ftcs_explicit_{tag}.txt")
    ferr_df  = os.path.join("data", f"error_dufort_{tag}.txt")
    ferr_imp = os.path.join("data", f"error_ftcs_implicit_{tag}.txt")
    ferr_cn  = os.path.join("data", f"error_cn_{tag}.txt")

    x = build_grid(imax, dx)
    L = (imax - 1) * dx
    alpha = F * (dx * dx) / dt

    Tex = np.zeros(imax + 1, dtype=float)
    Texp = np.zeros(imax + 1, dtype=float)
    Tdf = np.zeros(imax + 1, dtype=float)
    Timp = np.zeros(imax + 1, dtype=float)
    Tcn = np.zeros(imax + 1, dtype=float)

    with open(fexp, "w", encoding="utf-8") as out_exp, \
         open(fdf,  "w", encoding="utf-8") as out_df, \
         open(fimp, "w", encoding="utf-8") as out_imp, \
         open(fcn,  "w", encoding="utf-8") as out_cn, \
         open(fex,  "w", encoding="utf-8") as out_ex, \
         open(ferr_exp, "w", encoding="utf-8") as out_eexp, \
         open(ferr_df,  "w", encoding="utf-8") as out_edf, \
         open(ferr_imp, "w", encoding="utf-8") as out_eimp, \
         open(ferr_cn,  "w", encoding="utf-8") as out_ecn:

        for k in range(0, 5):
            t = 0.1 * k
            nmax = int(round(t / dt))

            FTCS_explicit(nmax, F, tb, t0, imax, Texp)
            DufortFrankel(nmax, F, tb, t0, imax, Tdf)
            FTCS_implicit(nmax, F, tb, t0, imax, Timp)
            CrankNicolson(nmax, F, tb, t0, imax, Tcn)

            exact_solution_profile(imax, x, t, alpha, L, tb, t0, 200, Tex)

            write_block_T(out_exp, t, x, Texp, imax)
            write_block_T(out_df,  t, x, Tdf,  imax)
            write_block_T(out_imp, t, x, Timp, imax)
            write_block_T(out_cn,  t, x, Tcn,  imax)
            write_block_T(out_ex,  t, x, Tex,  imax)

            write_block_error(out_eexp, t, x, Texp - Tex, imax)
            write_block_error(out_edf,  t, x, Tdf  - Tex, imax)
            write_block_error(out_eimp, t, x, Timp - Tex, imax)
            write_block_error(out_ecn,  t, x, Tcn  - Tex, imax)


# ============================================================================
# convergence: write one convergence table at t_target
# ============================================================================

def convergence_study(dx: float, imax: int, alpha: float,
                      t0: float, tb: float, t_target: float) -> str:
    dt_list = [0.10, 0.05, 0.02, 0.01, 0.005]
    outdata = os.path.join("data", f"convergence_t{make_t_tag(t_target)}.txt")

    x = build_grid(imax, dx)
    L = (imax - 1) * dx

    Tex = np.zeros(imax + 1, dtype=float)
    Texp = np.zeros(imax + 1, dtype=float)
    Timp = np.zeros(imax + 1, dtype=float)
    Tdf = np.zeros(imax + 1, dtype=float)
    Tcn = np.zeros(imax + 1, dtype=float)

    with open(outdata, "w", encoding="utf-8") as f:
        f.write(f"# convergence at t = {t_target:.2f} hr\n")
        f.write("# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\tLinf_exp\tLinf_imp\tLinf_df\tLinf_cn\n")

        for dt in dt_list:
            F = alpha * dt / (dx * dx)
            nmax = int(round(t_target / dt))

            FTCS_explicit(nmax, F, tb, t0, imax, Texp)
            FTCS_implicit(nmax, F, tb, t0, imax, Timp)
            DufortFrankel(nmax, F, tb, t0, imax, Tdf)
            CrankNicolson(nmax, F, tb, t0, imax, Tcn)

            exact_solution_profile(imax, x, t_target, alpha, L, tb, t0, 200, Tex)

            L2_exp = err_L2(imax, Texp, Tex)
            L2_imp = err_L2(imax, Timp, Tex)
            L2_df  = err_L2(imax, Tdf,  Tex)
            L2_cn  = err_L2(imax, Tcn,  Tex)

            Li_exp = err_Linf(imax, Texp, Tex)
            Li_imp = err_Linf(imax, Timp, Tex)
            Li_df  = err_Linf(imax, Tdf,  Tex)
            Li_cn  = err_Linf(imax, Tcn,  Tex)

            f.write(
                f"{dt:.6f}\t{L2_exp:.10e}\t{L2_imp:.10e}\t{L2_df:.10e}\t{L2_cn:.10e}\t"
                f"{Li_exp:.10e}\t{Li_imp:.10e}\t{Li_df:.10e}\t{Li_cn:.10e}\n"
            )

    return outdata






