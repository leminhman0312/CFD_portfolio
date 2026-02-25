module io_dirs
  implicit none
  private
  public:: ensure_dir, ensure_project_dirs, ensure_animation_dirs

contains
  subroutine ensure_dir(name)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: cmd
    integer :: exitstat
    cmd = 'mkdir -p ' // trim(name)
    call execute_command_line(cmd, exitstat=exitstat)
  end subroutine ensure_dir

  subroutine ensure_project_dirs()
    call ensure_dir('data')
    call ensure_dir('plot')
    call ensure_dir('gnuplot_scripts')
  end subroutine ensure_project_dirs

  subroutine ensure_animation_dirs()
    call ensure_dir('animation')
    call ensure_dir('animation/dat')
    call ensure_dir('animation/png')
    call ensure_dir('animation/movie')
  end subroutine ensure_animation_dirs
end module io_dirs

