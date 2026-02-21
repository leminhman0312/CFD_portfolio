module helpers
  use headers, only: dp
  implicit none
contains

  subroutine make_directory(dirname)
    implicit none
    character(len=*), intent(in) :: dirname
    integer :: exitstat

    call execute_command_line("mkdir -p " // trim(dirname), exitstat=exitstat)

    if (exitstat /= 0) then
      print *, "mkdir failed for: ", trim(dirname)
      stop 1
    end if
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
    do i = 1, imax - 1
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

end module helpers
