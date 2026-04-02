module io_field
  use kinds, only: real64
  implicit none
  private
  public :: write_field_xyz

contains

  subroutine write_field_xyz(filename, u, deltax, deltay)
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: u(:,:)
    real(real64), intent(in) :: deltax, deltay

    integer :: imax2, jmax2, i, j, unit
    real(real64) :: x, y

    imax2 = size(u,1)
    jmax2 = size(u,2)

    open(newunit=unit, file=trim(filename), status='replace', action='write')

    do j = 1, jmax2
      y = real(j-1, real64) * deltay
      do i = 1, imax2
        x = real(i-1, real64) * deltax
        write(unit,'(F12.8,1X,F12.8,1X,F16.10)') x, y, u(i,j)
      end do
      write(unit,*)
    end do

    close(unit)
  end subroutine write_field_xyz

end module io_field
