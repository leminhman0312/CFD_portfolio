module tridiag
  use kinds, only: real64
  implicit none
  private
  public :: thomasTriDiagonal

contains

  subroutine thomasTriDiagonal(n, a, b, c, d, u)
    integer, intent(in) :: n
    real(real64), intent(in) :: a(:), b(:), d(:)
    real(real64), intent(inout) :: c(:)
    real(real64), intent(inout) :: u(:)

    real(real64), allocatable :: dprime(:), cprime(:)
    integer :: i

    allocate(dprime(n), cprime(n))

    dprime(1) = d(1)
    cprime(1) = c(1)

    do i = 2, n
      dprime(i) = d(i) - (b(i) * a(i-1)) / dprime(i-1)
      cprime(i) = c(i) - (cprime(i-1) * b(i)) / dprime(i-1)
    end do

    u(n) = cprime(n) / dprime(n)

    do i = n-1, 1, -1
      u(i) = (cprime(i) - a(i) * u(i+1)) / dprime(i)
    end do

    deallocate(dprime, cprime)
  end subroutine thomasTriDiagonal

end module tridiag
