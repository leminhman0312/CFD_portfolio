module methods
  use headers, only: dp
  implicit none
contains

  ! norms
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

  ! thomas solver
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

    do i = N - 1, 1, -1
      u(i) = (cprime(i) - a(i) * u(i+1)) / dprime(i)
    end do

    deallocate(dprime)
    deallocate(cprime)
  end subroutine thomas_solver

  ! exact solution
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

        ! For even n: (1 - (-1)^n) = 0, term is exactly zero
        if (mod(n, 2) == 0) cycle

        ! For odd n: (1 - (-1)^n) = 2
        coeff = (2.0_dp * A / (real(n, dp) * pi)) * 2.0_dp
        k = real(n, dp) * pi / L

        s = s + coeff * sin(k * xi) * exp(-alpha * k * k * t)
      end do

      Tex(i) = tb + s
    end do

    Tex(1)    = tb
    Tex(imax) = tb
  end subroutine exact_solution

  ! schemes
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
    do i = 2, imax - 1
      un(i) = t0
    end do

    do n = 1, nmax
      u = un
      do i = 2, imax - 1
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
    do i = 2, imax - 1
      un(i) = t0
    end do

    d(1) = 1.0_dp
    d(imax) = 1.0_dp
    a(1) = 0.0_dp
    b(1) = 0.0_dp
    a(imax) = 0.0_dp
    b(imax) = 0.0_dp

    do i = 2, imax - 1
      d(i) = 1.0_dp + 2.0_dp*F
      a(i) = -F
      b(i) = -F
    end do

    do step = 1, nmax
      c(1) = tb
      c(imax) = tb
      do i = 2, imax - 1
        c(i) = un(i)
      end do

      call thomas_solver(imax, a, b, c, d, un)

      un(1) = tb
      un(imax) = tb
    end do

    Tout = un

    deallocate(a, b, c, d, un)
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

    allocate(a(1:imax), b(1:imax), c(1:imax), d(1:imax))
    allocate(u0(1:imax), uhalf(1:imax))

    u0(1) = tb
    u0(imax) = tb
    do i = 2, imax - 1
      u0(i) = t0
    end do

    d(1) = 1.0_dp
    d(imax) = 1.0_dp
    a(1) = 0.0_dp
    b(1) = 0.0_dp
    a(imax) = 0.0_dp
    b(imax) = 0.0_dp

    do i = 2, imax - 1
      d(i) = 1.0_dp + 2.0_dp*dnc
      a(i) = -dnc
      b(i) = -dnc
    end do

    do step = 1, nmax

      uhalf(1) = tb
      uhalf(imax) = tb
      do i = 2, imax - 1
        uhalf(i) = u0(i) + dnc * (u0(i+1) - 2.0_dp*u0(i) + u0(i-1))
      end do

      c(1) = tb
      c(imax) = tb
      do i = 2, imax - 1
        c(i) = uhalf(i)
      end do

      call thomas_solver(imax, a, b, c, d, u0)

      u0(1) = tb
      u0(imax) = tb
    end do

    Tout = u0

    deallocate(a, b, c, d, u0, uhalf)
  end subroutine scheme_crank_nicolson

  subroutine implicit_interior_one_step(F, tb, t0, imax, u_full)
    implicit none
    integer, intent(in) :: imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: u_full(1:imax)

    integer :: N, j
    real(dp), allocatable :: a(:), b(:), c(:), d(:), u_int(:)

    N = imax - 2

    ! Safety guard for tiny grids
    if (N < 1) then
      u_full(1) = tb
      u_full(imax) = tb
      return
    end if

    allocate(a(1:N), b(1:N), c(1:N), d(1:N), u_int(1:N))

    ! Make compiler and runtime happy: everything initialized
    a = 0.0_dp
    b = 0.0_dp
    c = 0.0_dp
    d = 0.0_dp
    u_int = 0.0_dp

    u_full(1) = tb
    u_full(imax) = tb
    do j = 1, N
      u_full(j+1) = t0
    end do

    do j = 1, N
      d(j) = 1.0_dp + 2.0_dp * F
      a(j) = -F
      b(j) = -F
    end do

    b(1) = 0.0_dp
    a(N) = 0.0_dp

    do j = 1, N
      c(j) = u_full(j+1)
    end do

    c(1) = c(1) + F * tb
    c(N) = c(N) + F * tb

    call thomas_solver(N, a, b, c, d, u_int)

    u_full(1) = tb
    u_full(imax) = tb
    do j = 1, N
      u_full(j+1) = u_int(j)
    end do

    deallocate(a, b, c, d, u_int)
  end subroutine implicit_interior_one_step

  subroutine scheme_dufort_frankel(nmax, F, tb, t0, imax, Tout)
    implicit none
    integer, intent(in) :: nmax, imax
    real(dp), intent(in) :: F, tb, t0
    real(dp), intent(out) :: Tout(1:imax)

    real(dp) :: dval
    real(dp), allocatable :: un_m1(:), un(:), un_p1(:), u1(:)
    integer :: i, step

    dval = 2.0_dp * F

    allocate(un_m1(1:imax), un(1:imax), un_p1(1:imax), u1(1:imax))

    un_m1(1) = tb
    un_m1(imax) = tb
    do i = 2, imax - 1
      un_m1(i) = t0
    end do

    call implicit_interior_one_step(F, tb, t0, imax, u1)
    un = u1

    if (nmax <= 0) then
      Tout = un_m1
      deallocate(un_m1, un, un_p1, u1)
      return
    end if

    do step = 1, nmax - 1
      do i = 2, imax - 1
        un_p1(i) = ((1.0_dp - dval) * un_m1(i) + dval * (un(i+1) + un(i-1))) / (1.0_dp + dval)
      end do
      un_p1(1) = tb
      un_p1(imax) = tb

      un_m1 = un
      un    = un_p1
    end do

    Tout = un

    deallocate(un_m1, un, un_p1, u1)
  end subroutine scheme_dufort_frankel

  ! io blocks
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

end module methods
