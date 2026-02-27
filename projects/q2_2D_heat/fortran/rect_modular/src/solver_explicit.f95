module solvers_explicit
  use kinds, only: real64
  implicit none
  private
  public :: FTCS_Explicit

contains

  subroutine FTCS_Explicit(u0, t_end, deltax, deltay, dt, alpha, t1, t2, t3, t4, u)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, dt, alpha
    real(real64), intent(in) :: t1, t2, t3, t4
    real(real64), intent(out) :: u(size(u0,1), size(u0,2))

    integer :: imax2, jmax2
    integer :: i, j, n, nfull, nsteps
    real(real64) :: dt_last, dt_n, fx, fy
    real(real64), allocatable :: u_new(:,:)

    imax2 = size(u0,1)
    jmax2 = size(u0,2)

    allocate(u_new(imax2, jmax2))
    u = u0

    nfull = int(floor(t_end / dt))
    dt_last = t_end - real(nfull, real64) * dt
    if (dt_last > 0.0_real64) then
      nsteps = nfull + 1
    else
      nsteps = nfull
    end if

    do n = 1, nsteps

      if (n == nsteps .and. dt_last > 0.0_real64) then
        dt_n = dt_last
      else
        dt_n = dt
      end if

      fx = alpha * dt_n / (deltax*deltax)
      fy = alpha * dt_n / (deltay*deltay)

      u_new = u

      do i = 2, imax2-1
        do j = 2, jmax2-1
          u_new(i,j) = (1.0_real64 - 2.0_real64*fx - 2.0_real64*fy) * u(i,j) &
                     + fx * (u(i+1,j) + u(i-1,j)) &
                     + fy * (u(i,j+1) + u(i,j-1))
        end do
      end do

      do j = 1, jmax2
        u_new(1,j)     = t2
        u_new(imax2,j) = t4
      end do
      do i = 1, imax2
        u_new(i,1)     = t1
        u_new(i,jmax2) = t3
      end do

      u = u_new
    end do

    deallocate(u_new)
  end subroutine FTCS_Explicit

end module solvers_explicit
