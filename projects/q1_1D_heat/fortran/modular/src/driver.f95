module driver
  use headers, only: dp
  use helpers, only: make_dt_tag, make_t_tag, build_grid, str_f06, str_f02, str_f01, str_i0
  use methods, only: exact_solution, &
                     scheme_ftcs_explicit, scheme_ftcs_implicit, scheme_dufort_frankel, scheme_crank_nicolson, &
                     write_block_T, write_block_error, &
                     error_L2, error_Linf
  implicit none
contains

  subroutine sim(dx, dt, imax, t0, tb, F)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: dx, dt, t0, tb, F

    character(len=3) :: dtag
    character(len=256) :: fexp, fimp, fdf, fcn, fex
    character(len=256) :: feexp, feimp, fedf, fecn

    integer :: uexp, uimp, udf, ucn, uex
    integer :: ueexp, ueimp, uedf, uecn

    real(dp), allocatable :: x(:)
    real(dp), allocatable :: Tex(:), Texp(:), Timp(:), Tdf(:), Tcn(:)
    real(dp) :: L, alpha_loc, t
    integer :: k, nmax

    call make_dt_tag(dt, dtag)

    fexp  = "data/ftcs_explicit_" // dtag // ".txt"
    fimp  = "data/ftcs_implicit_" // dtag // ".txt"
    fdf   = "data/dufort_"        // dtag // ".txt"
    fcn   = "data/cn_"            // dtag // ".txt"
    fex   = "data/exact_"         // dtag // ".txt"

    feexp = "data/error_ftcs_explicit_" // dtag // ".txt"
    feimp = "data/error_ftcs_implicit_" // dtag // ".txt"
    fedf  = "data/error_dufort_"        // dtag // ".txt"
    fecn  = "data/error_cn_"            // dtag // ".txt"

    allocate(x(1:imax), Tex(1:imax), Texp(1:imax), Timp(1:imax), Tdf(1:imax), Tcn(1:imax))

    call build_grid(imax, dx, x)

    L = real(imax-1, dp) * dx
    alpha_loc = F * (dx*dx) / dt

    uexp  = 10
    uimp  = 11
    udf   = 12
    ucn   = 13
    uex   = 14

    ueexp = 20
    ueimp = 21
    uedf  = 22
    uecn  = 23

    open(uexp,  file=trim(fexp),  status="replace", action="write")
    open(uimp,  file=trim(fimp),  status="replace", action="write")
    open(udf,   file=trim(fdf),   status="replace", action="write")
    open(ucn,   file=trim(fcn),   status="replace", action="write")
    open(uex,   file=trim(fex),   status="replace", action="write")

    open(ueexp, file=trim(feexp), status="replace", action="write")
    open(ueimp, file=trim(feimp), status="replace", action="write")
    open(uedf,  file=trim(fedf),  status="replace", action="write")
    open(uecn,  file=trim(fecn),  status="replace", action="write")

    do k = 0, 4
      t = 0.1_dp * real(k, dp)
      nmax = nint(t / dt)

      call scheme_ftcs_explicit(nmax, F, tb, t0, imax, Texp)
      call scheme_ftcs_implicit(nmax, F, tb, t0, imax, Timp)
      call scheme_dufort_frankel(nmax, F, tb, t0, imax, Tdf)
      call scheme_crank_nicolson(nmax, F, tb, t0, imax, Tcn)

      call exact_solution(imax, x, t, alpha_loc, L, tb, t0, 200, Tex)

      call write_block_T(uexp, t, imax, x, Texp)
      call write_block_T(uimp, t, imax, x, Timp)
      call write_block_T(udf,  t, imax, x, Tdf)
      call write_block_T(ucn,  t, imax, x, Tcn)
      call write_block_T(uex,  t, imax, x, Tex)

      call write_block_error(ueexp, t, imax, x, Texp, Tex)
      call write_block_error(ueimp, t, imax, x, Timp, Tex)
      call write_block_error(uedf,  t, imax, x, Tdf,  Tex)
      call write_block_error(uecn,  t, imax, x, Tcn,  Tex)
    end do

    close(uexp);  close(uimp);  close(udf);  close(ucn);  close(uex)
    close(ueexp); close(ueimp); close(uedf); close(uecn)

    deallocate(x, Tex, Texp, Timp, Tdf, Tcn)
  end subroutine sim

  subroutine convergence_study(dx, imax, alpha, t0, tb, t_target)
    implicit none
    real(dp), intent(in) :: dx, alpha, t0, tb, t_target
    integer, intent(in) :: imax

    integer, parameter :: ndt = 5
    real(dp) :: dt_list(ndt)
    real(dp) :: dt, F, L
    integer :: j, nmax, uout
    character(len=3) :: ttag
    character(len=256) :: outdata

    real(dp), allocatable :: x(:)
    real(dp), allocatable :: Tex(:), Texp(:), Timp(:), Tdf(:), Tcn(:)

    real(dp) :: L2e, L2i, L2d, L2c
    real(dp) :: Lie, Lii, Lid, Lic

    dt_list = [0.10_dp, 0.05_dp, 0.02_dp, 0.01_dp, 0.005_dp]

    call make_t_tag(t_target, ttag)
    outdata = "data/convergence_t" // ttag // ".txt"

    allocate(x(1:imax), Tex(1:imax), Texp(1:imax), Timp(1:imax), Tdf(1:imax), Tcn(1:imax))

    call build_grid(imax, dx, x)
    L = real(imax-1, dp) * dx

    uout = 80
    open(uout, file=trim(outdata), status="replace", action="write")
    write(uout, "(A,F6.2,A)") "# convergence at t = ", t_target, " hr"
    write(uout, "(A)") "# dt<TAB>L2_exp<TAB>L2_imp<TAB>L2_df<TAB>L2_cn<TAB>Linf_exp<TAB>Linf_imp<TAB>Linf_df<TAB>Linf_cn"

    do j = 1, ndt
      dt = dt_list(j)
      F  = alpha * dt / (dx*dx)
      nmax = nint(t_target / dt)

      call scheme_ftcs_explicit(nmax, F, tb, t0, imax, Texp)
      call scheme_ftcs_implicit(nmax, F, tb, t0, imax, Timp)
      call scheme_dufort_frankel(nmax, F, tb, t0, imax, Tdf)
      call scheme_crank_nicolson(nmax, F, tb, t0, imax, Tcn)

      call exact_solution(imax, x, t_target, alpha, L, tb, t0, 200, Tex)

      L2e = error_L2(imax, Texp, Tex)
      L2i = error_L2(imax, Timp, Tex)
      L2d = error_L2(imax, Tdf,  Tex)
      L2c = error_L2(imax, Tcn,  Tex)

      Lie = error_Linf(imax, Texp, Tex)
      Lii = error_Linf(imax, Timp, Tex)
      Lid = error_Linf(imax, Tdf,  Tex)
      Lic = error_Linf(imax, Tcn,  Tex)

      write(uout, "(F10.6,A,8(ES16.8,A))") dt, achar(9), &
        L2e, achar(9), L2i, achar(9), L2d, achar(9), L2c, achar(9), &
        Lie, achar(9), Lii, achar(9), Lid, achar(9), Lic, achar(9)
    end do

    close(uout)
    deallocate(x, Tex, Texp, Timp, Tdf, Tcn)
  end subroutine convergence_study

  subroutine plot(dt, t_target)
    implicit none
    real(dp), intent(in) :: dt, t_target

    character(len=3) :: tag, ttag
    character(len=512) :: cmd
    character(len=256) :: infile
    integer :: idx

    call make_dt_tag(dt, tag)
    call make_t_tag(t_target, ttag)

    cmd = "gnuplot -e ""scheme='ftcs_explicit'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/ftcs_explicit_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call execute_command_line(trim(cmd))

    cmd = "gnuplot -e ""scheme='dufort'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/dufort_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call execute_command_line(trim(cmd))

    cmd = "gnuplot -e ""scheme='ftcs_implicit'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/ftcs_implicit_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call execute_command_line(trim(cmd))

    cmd = "gnuplot -e ""scheme='cn'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/cn_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call execute_command_line(trim(cmd))

    idx = nint(t_target / 0.1_dp)

    cmd = "gnuplot -e ""tag='" // tag // "'; dt=" // trim(str_f06(dt)) // &
          "; idx=" // trim(str_i0(idx)) // &
          "; tlabel='" // trim(str_f01(t_target)) // &
          "'; outpng='plot/error_schemes_" // tag // "_t" // trim(str_f01(t_target)) // ".png'"" " // &
          "gnuplot_scripts/compare_error_schemes.gp"

    call execute_command_line(trim(cmd))
    
    infile = "data/convergence_t" // ttag // ".txt"
    cmd = "gnuplot -e ""infile='" // trim(infile) // &
          "'; tlabel='" // trim(str_f02(t_target)) // &
          "'; outpng='plot/convergence_t" // trim(str_f02(t_target)) // ".png'"" " // &
          "gnuplot_scripts/convergence.gp"
    call execute_command_line(trim(cmd))
  end subroutine plot

end module driver
