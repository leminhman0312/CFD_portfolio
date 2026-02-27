module field_init
  use kinds, only: real64
  implicit none
  private
  public :: initializeField
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

end module field_init












