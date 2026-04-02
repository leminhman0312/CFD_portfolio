module field_init
  use kinds, only: real64
  implicit none
  private
  public :: initializeField, initialize_point_source
contains

  subroutine initializeField(imax, jmax, t0, t1, t2, t3, t4, u)
    integer, intent(in) :: imax, jmax
    real(real64), intent(in) :: t0, t1, t2, t3, t4
    real(real64), intent(out) :: u(imax, jmax)
    integer :: i, j

    u = t0

    do j = 1, jmax
      u(1, j)    = t2
      u(imax, j) = t4
    end do
    do i = 1, imax
      u(i, 1)    = t1
      u(i, jmax) = t3
    end do
  end subroutine initializeField

  subroutine initialize_point_source(imax, jmax, deltax, deltay, r_source, t_source, u)
    integer, intent(in) :: imax, jmax
    real(real64), intent(in) :: deltax, deltay, r_source, t_source
    real(real64), intent(out) :: u(imax, jmax)

    u = 0.0_real64
    call clamp_circle(u, deltax, deltay, r_source, t_source)
    call apply_bc_zero(u)
  end subroutine initialize_point_source

  subroutine clamp_circle(u, deltax, deltay, r, tval)
    real(real64), intent(inout) :: u(:,:)
    real(real64), intent(in) :: deltax, deltay, r, tval
    integer :: im, jm, i, j
    real(real64) :: x0, y0, x, y, dx, dy, rr2

    im = size(u,1)
    jm = size(u,2)

    x0 = 0.5_real64 * real(im - 1, real64) * deltax
    y0 = 0.5_real64 * real(jm - 1, real64) * deltay
    rr2 = r*r

    do j = 1, jm
      y  = real(j - 1, real64) * deltay
      dy = y - y0
      do i = 1, im
        x  = real(i - 1, real64) * deltax
        dx = x - x0
        if (dx*dx + dy*dy <= rr2) u(i,j) = tval
      end do
    end do
  end subroutine clamp_circle

  subroutine apply_bc_zero(u)
    real(real64), intent(inout) :: u(:,:)
    integer :: im, jm, i, j

    im = size(u,1)
    jm = size(u,2)

    do j = 1, jm
      u(1,j)  = 0.0_real64
      u(im,j) = 0.0_real64
    end do
    do i = 1, im
      u(i,1)  = 0.0_real64
      u(i,jm) = 0.0_real64
    end do
  end subroutine apply_bc_zero

end module field_init












