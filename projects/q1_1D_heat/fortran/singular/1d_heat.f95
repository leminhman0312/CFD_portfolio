! ============================================================================
! heat_1d.f95
!
! Goal
!   A very easy to follow single Fortran file you can learn from
!
! What it does (matches the clean Cpp behavior)
!   1) Runs 1D heat equation with 4 schemes
!      FTCS explicit
!      FTCS implicit
!      Dufort Frankel
!      Crank Nicolson
!   2) Writes data files with 5 time blocks t = 0.0 0.1 0.2 0.3 0.4
!      Also writes error files (scheme minus exact) for the same blocks
!   3) Runs a convergence study and writes one table file
!   4) Calls your existing gnuplot scripts and writes png into plot
!
! Compile
!   gfortran -O2 heat_1d.f95 -o heat1d
!
! Run
!   ./heat1d
! ============================================================================

program heat_1d
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)

  integer  :: imax
  real(dp) :: dx, dt, alpha, t0, tb, Fourier_num
  real(dp) :: t_target

  imax  = 21
  dx    = 0.05_dp
  dt    = 0.01_dp
  alpha = 0.1_dp
  t0    = 100.0_dp
  tb    = 300.0_dp

  Fourier_num = alpha * dt / (dx*dx)

  t_target = 0.4_dp

  call make_directory("data")
  call make_directory("plot")

  call sim(dx, dt, imax, t0, tb, Fourier_num)

  call convergence_study(dx, imax, alpha, t0, tb, t_target)

  call plot(dt, t_target)

  print *, "Done"

contains

  ! ==========================================================================
  ! Section A
  ! Small helpers
  ! ==========================================================================

  subroutine make_directory(dirname)
    implicit none
    character(len=*), intent(in) :: dirname
    call system("mkdir -p " // trim(dirname))
  end subroutine make_directory

  subroutine make_dt_tag(dt, tag)
    implicit none
    real(dp), intent(in) :: dt
    character(len=3), intent(out) :: tag
    integer :: code
    code = nint(dt * 100.0_dp)
    write(tag, "(I3.3)") code
  end subroutine make_dt_tag

  subroutine make_t_tag(t, tag)
    implicit none
    real(dp), intent(in) :: t
    character(len=3), intent(out) :: tag
    integer :: code
    code = nint(t * 100.0_dp)
    write(tag, "(I3.3)") code
  end subroutine make_t_tag

  subroutine build_grid(imax, dx, x)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: x(1:imax)
    integer :: i
    x(1) = 0.0_dp
    do i = 1, imax-1
      x(i+1) = x(i) + dx
    end do
  end subroutine build_grid

  function str_f06(x) result(s)
    implicit none
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s, "(F0.6)") x
    s = adjustl(s)
  end function str_f06

  function str_f02(x) result(s)
    implicit none
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s, "(F0.2)") x
    s = adjustl(s)
  end function str_f02

  function str_f01(x) result(s)
    implicit none
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s, "(F0.1)") x
    s = adjustl(s)
  end function str_f01

  function str_i0(i) result(s)
    implicit none
    integer, intent(in) :: i
    character(len=32) :: s
    write(s, "(I0)") i
    s = adjustl(s)
  end function str_i0

  ! ==========================================================================
  ! Section B
  ! Error norms (used in convergence table)
  ! ==========================================================================

  real(dp) function error_L2(imax, a, b)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: a(1:imax), b(1:imax)
    integer :: i
    real(dp) :: s, e
    s = 0.0_dp
    do i = 1, imax
      e = a(i) - b(i)
      s = s + e*e
    end do
    error_L2 = sqrt(s / real(imax, dp))
  end function error_L2

  real(dp) function error_Linf(imax, a, b)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: a(1:imax), b(1:imax)
    integer :: i
    real(dp) :: m, e
    m = 0.0_dp
    do i = 1, imax
      e = abs(a(i) - b(i))
      if (e > m) m = e
    end do
    error_Linf = m
  end function error_Linf

  ! ==========================================================================
  ! Section C
  ! Tridiagonal solver
  !
  ! Solves for u(1:N) given a,b,c,d arrays indexed 1:N.
  ! This matches the Thomas algorithm used in the clean Cpp.
  ! ==========================================================================

  subroutine thomas_solver(N, a, b, c, d, u)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in)  :: a(1:N), b(1:N), c(1:N), d(1:N)
    real(dp), intent(out) :: u(1:N)

    real(dp), allocatable :: dprime(:), cprime(:)
    integer :: i

    allocate(dprime(1:N))
    allocate(cprime(1:N))

    dprime(1) = d(1)
    cprime(1) = c(1)

    do i = 2, N
      dprime(i) = d(i) - (b(i) * a(i-1)) / dprime(i-1)
      cprime(i) = c(i) - (cprime(i-1) * b(i)) / dprime(i-1)
    end do

    u(N) = cprime(N) / dprime(N)

    do i = N-1, 1, -1
      u(i) = (cprime(i) - a(i) * u(i+1)) / dprime(i)
    end do

    deallocate(dprime)
    deallocate(cprime)
  end subroutine thomas_solver

  ! ==========================================================================
  ! Section D
  ! Exact solution
  ! ==========================================================================

  subroutine exact_solution(imax, x, t, alpha, L, tb, t0, nterms, Tex)
    implicit none
    integer, intent(in) :: imax, nterms
    real(dp), intent(in) :: x(1:imax), t, alpha, L, tb, t0
    real(dp), intent(out) :: Tex(1:imax)

    integer :: i, n
    real(dp) :: pi, A, s, xi, coeff, k

    pi = 4.0_dp * atan(1.0_dp)
    A  = t0 - tb

    do i = 1, imax
      xi = x(i)
      s  = 0.0_dp
      do n = 1, nterms
        coeff = (2.0_dp * A / (real(n,dp) * pi)) * (1.0_dp - (-1.0_dp)**n)
        if (coeff == 0.0_dp) cycle
        k = real(n,dp) * pi / L
        s = s + coeff * sin(k*xi) * exp(-alpha * k*k * t)
      end do
      Tex(i) = tb + s
    end do

    Tex(1)    = tb
    Tex(imax) = tb
  end subroutine exact_solution

  ! ==========================================================================
  ! Section E
  ! Numerical schemes
  ! Each scheme fills Tout(1:imax)
  ! ==========================================================================

  subroutine scheme_ftcs_explicit(nmax, F, tb, t0, imax, Tout)
    implicit none
    integer, intent(in) :: nmax, imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: Tout(1:imax)

    real(dp), allocatable :: u(:), un(:)
    integer :: i, n

    allocate(u(1:imax))
    allocate(un(1:imax))

    un(1) = tb
    un(imax) = tb
    do i = 2, imax-1
      un(i) = t0
    end do

    do n = 1, nmax
      u = un
      do i = 2, imax-1
        un(i) = u(i) + F * (u(i+1) - 2.0_dp*u(i) + u(i-1))
      end do
      un(1) = tb
      un(imax) = tb
    end do

    Tout = un

    deallocate(u)
    deallocate(un)
  end subroutine scheme_ftcs_explicit

  subroutine scheme_ftcs_implicit(nmax, F, tb, t0, imax, Tout)
    implicit none
    integer, intent(in) :: nmax, imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: Tout(1:imax)

    real(dp), allocatable :: a(:), b(:), c(:), d(:), un(:)
    integer :: i, step

    allocate(a(1:imax))
    allocate(b(1:imax))
    allocate(c(1:imax))
    allocate(d(1:imax))
    allocate(un(1:imax))

    un(1) = tb
    un(imax) = tb
    do i = 2, imax-1
      un(i) = t0
    end do

    d(1) = 1.0_dp
    d(imax) = 1.0_dp
    a(1) = 0.0_dp
    b(1) = 0.0_dp
    a(imax) = 0.0_dp
    b(imax) = 0.0_dp

    do i = 2, imax-1
      d(i) = 1.0_dp + 2.0_dp*F
      a(i) = -F
      b(i) = -F
    end do

    do step = 1, nmax
      c(1) = tb
      c(imax) = tb
      do i = 2, imax-1
        c(i) = un(i)
      end do

      call thomas_solver(imax, a, b, c, d, un)

      un(1) = tb
      un(imax) = tb
    end do

    Tout = un

    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(d)
    deallocate(un)
  end subroutine scheme_ftcs_implicit

  subroutine scheme_crank_nicolson(nmax, F, tb, t0, imax, Tout)
    implicit none
    integer, intent(in) :: nmax, imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: Tout(1:imax)

    real(dp) :: dnc
    real(dp), allocatable :: a(:), b(:), c(:), d(:), u0(:), uhalf(:)
    integer :: i, step

    dnc = F / 2.0_dp

    allocate(a(1:imax))
    allocate(b(1:imax))
    allocate(c(1:imax))
    allocate(d(1:imax))
    allocate(u0(1:imax))
    allocate(uhalf(1:imax))

    u0(1) = tb
    u0(imax) = tb
    do i = 2, imax-1
      u0(i) = t0
    end do

    d(1) = 1.0_dp
    d(imax) = 1.0_dp
    a(1) = 0.0_dp
    b(1) = 0.0_dp
    a(imax) = 0.0_dp
    b(imax) = 0.0_dp

    do i = 2, imax-1
      d(i) = 1.0_dp + 2.0_dp*dnc
      a(i) = -dnc
      b(i) = -dnc
    end do

    do step = 1, nmax

      uhalf(1) = tb
      uhalf(imax) = tb
      do i = 2, imax-1
        uhalf(i) = u0(i) + dnc * (u0(i+1) - 2.0_dp*u0(i) + u0(i-1))
      end do

      c(1) = tb
      c(imax) = tb
      do i = 2, imax-1
        c(i) = uhalf(i)
      end do

      call thomas_solver(imax, a, b, c, d, u0)

      u0(1) = tb
      u0(imax) = tb
    end do

    Tout = u0

    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(d)
    deallocate(u0)
    deallocate(uhalf)
  end subroutine scheme_crank_nicolson

  subroutine scheme_dufort_frankel(nmax, F, tb, t0, imax, Tout)
    implicit none
    integer, intent(in) :: nmax, imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: Tout(1:imax)

    real(dp) :: dval
    real(dp), allocatable :: un_m1(:), un(:), un_p1(:), u1(:)
    integer :: i, step

    dval = 2.0_dp * F

    allocate(un_m1(1:imax))
    allocate(un(1:imax))
    allocate(un_p1(1:imax))
    allocate(u1(1:imax))

    un_m1(1) = tb
    un_m1(imax) = tb
    do i = 2, imax-1
      un_m1(i) = t0
    end do

    call implicit_interior_one_step(F, tb, t0, imax, u1)
    un = u1

    if (nmax <= 0) then
      Tout = un_m1
      deallocate(un_m1, un, un_p1, u1)
      return
    end if

    do step = 1, nmax-1
      do i = 2, imax-1
        un_p1(i) = ((1.0_dp - dval) * un_m1(i) + dval * (un(i+1) + un(i-1))) / (1.0_dp + dval)
      end do
      un_p1(1) = tb
      un_p1(imax) = tb

      un_m1 = un
      un    = un_p1
    end do

    Tout = un

    deallocate(un_m1)
    deallocate(un)
    deallocate(un_p1)
    deallocate(u1)
  end subroutine scheme_dufort_frankel

  subroutine implicit_interior_one_step(F, tb, t0, imax, u_full)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: u_full(1:imax)

    integer :: N, j
    real(dp), allocatable :: a(:), b(:), c(:), d(:), u_int(:)

    N = imax - 2

    allocate(a(1:N))
    allocate(b(1:N))
    allocate(c(1:N))
    allocate(d(1:N))
    allocate(u_int(1:N))

    u_full(1) = tb
    u_full(imax) = tb
    do j = 1, N
      u_full(j+1) = t0
    end do

    do j = 1, N
      d(j) = 1.0_dp + 2.0_dp*F
      a(j) = -F
      b(j) = -F
    end do

    b(1) = 0.0_dp
    a(N) = 0.0_dp

    do j = 1, N
      c(j) = u_full(j+1)
    end do
    c(1) = c(1) + F*tb
    c(N) = c(N) + F*tb

    call thomas_solver(N, a, b, c, d, u_int)

    u_full(1) = tb
    u_full(imax) = tb
    do j = 1, N
      u_full(j+1) = u_int(j)
    end do

    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(d)
    deallocate(u_int)
  end subroutine implicit_interior_one_step

  ! ==========================================================================
  ! Section F
  ! Writing simulation files
  ! ==========================================================================

  subroutine write_block_T(unit, t, imax, x, Tarr)
    implicit none
    integer, intent(in) :: unit, imax
    real(dp), intent(in) :: t
    real(dp), intent(in) :: x(1:imax), Tarr(1:imax)
    integer :: i

    write(unit, "(A,F5.2,A)") "# t = ", t, " hr"
    write(unit, "(A)") "# x<TAB>T"
    do i = 1, imax
      write(unit, "(F10.6,A,F10.6)") x(i), achar(9), Tarr(i)
    end do
    write(unit,*) ""
    write(unit,*) ""
  end subroutine write_block_T

  subroutine write_block_error(unit, t, imax, x, Tarr, Tex)
    implicit none
    integer, intent(in) :: unit, imax
    real(dp), intent(in) :: t
    real(dp), intent(in) :: x(1:imax), Tarr(1:imax), Tex(1:imax)
    integer :: i

    write(unit, "(A,F5.2,A)") "# t = ", t, " hr"
    write(unit, "(A)") "# x<TAB>err"
    do i = 1, imax
      write(unit, "(F10.6,A,F10.6)") x(i), achar(9), (Tarr(i) - Tex(i))
    end do
    write(unit,*) ""
    write(unit,*) ""
  end subroutine write_block_error

  ! ==========================================================================
  ! Section G
  ! sim
  ! Writes solution files and error files for one dt
  ! ==========================================================================

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

    allocate(x(1:imax))
    allocate(Tex(1:imax))
    allocate(Texp(1:imax))
    allocate(Timp(1:imax))
    allocate(Tdf(1:imax))
    allocate(Tcn(1:imax))

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

    close(uexp)
    close(uimp)
    close(udf)
    close(ucn)
    close(uex)

    close(ueexp)
    close(ueimp)
    close(uedf)
    close(uecn)

    deallocate(x)
    deallocate(Tex)
    deallocate(Texp)
    deallocate(Timp)
    deallocate(Tdf)
    deallocate(Tcn)
  end subroutine sim

  ! ==========================================================================
  ! Section H
  ! convergence_study
  ! Writes data/convergence_tTTT.txt where TTT = t_target*100
  ! ==========================================================================

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

    dt_list(1) = 0.10_dp
    dt_list(2) = 0.05_dp
    dt_list(3) = 0.02_dp
    dt_list(4) = 0.01_dp
    dt_list(5) = 0.005_dp

    call make_t_tag(t_target, ttag)
    outdata = "data/convergence_t" // ttag // ".txt"

    allocate(x(1:imax))
    allocate(Tex(1:imax))
    allocate(Texp(1:imax))
    allocate(Timp(1:imax))
    allocate(Tdf(1:imax))
    allocate(Tcn(1:imax))

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

    deallocate(x)
    deallocate(Tex)
    deallocate(Texp)
    deallocate(Timp)
    deallocate(Tdf)
    deallocate(Tcn)
  end subroutine convergence_study

  ! ==========================================================================
  ! Section I
  ! plot
  ! Calls your shared gnuplot scripts and writes png into plot
  ! ==========================================================================

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
    call system(trim(cmd))

    cmd = "gnuplot -e ""scheme='dufort'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/dufort_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call system(trim(cmd))

    cmd = "gnuplot -e ""scheme='ftcs_implicit'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/ftcs_implicit_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call system(trim(cmd))

    cmd = "gnuplot -e ""scheme='cn'; tag='" // tag // &
          "'; dt='" // trim(str_f06(dt)) // &
          "'; outpng='plot/cn_" // tag // ".png'"" gnuplot_scripts/plot_scheme.gp"
    call system(trim(cmd))

    idx = nint(t_target / 0.1_dp)

    cmd = "gnuplot -e ""tag='" // tag // "'; dt=" // trim(str_f06(dt)) // &
          "; idx=" // trim(str_i0(idx)) // &
          "; tlabel='" // trim(str_f01(t_target)) // &
          "'; outpng='plot/error_schemes_" // tag // "_t" // trim(str_f01(t_target)) // ".png'"" " // &
          "gnuplot_scripts/compare_error_schemes.gp"
    call system(trim(cmd))

    infile = "data/convergence_t" // ttag // ".txt"
    cmd = "gnuplot -e ""infile='" // trim(infile) // &
          "'; tlabel='" // trim(str_f02(t_target)) // &
          "'; outpng='plot/convergence_t" // trim(str_f02(t_target)) // ".png'"" " // &
          "gnuplot_scripts/convergence.gp"
    call system(trim(cmd))
  end subroutine plot

end program heat_1d
