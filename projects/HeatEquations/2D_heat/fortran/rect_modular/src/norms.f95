module norms
  use kinds, only: real64
  implicit none
  private
  public :: error_L2, error_Linf

contains
  
  real(real64) function error_L2(u, uref) result(val)
    real(real64), intent(in) :: u(:,:), uref(:,:)
    integer :: imax, jmax, i, j, n
    real(real64) :: sum, e

    imax = size(u,1)
    jmax = size(u,2)

    sum = 0.0_real64
    n = 0

    do i = 2, imax-1
      do j = 2, jmax-1
        e = u(i,j) - uref(i,j)
        sum = sum + e*e
        n = n + 1
      end do
    end do

    if (n > 0) then
      val = sqrt(sum / real(n, real64))
    else
      val = 0.0_real64
    end if
  end function error_L2

  real(real64) function error_Linf(u, uref) result(val)
    real(real64), intent(in) :: u(:,:), uref(:,:)
    integer :: imax, jmax, i, j
    real(real64) :: e

    imax = size(u,1)
    jmax = size(u,2)

    val = 0.0_real64
    do i = 2, imax-1
      do j = 2, jmax-1
        e = abs(u(i,j) - uref(i,j))
        if (e > val) val = e
      end do
    end do
  end function error_Linf

end module norms
