module animation
  use kinds,       only: real64
  use solvers_adi, only: FTCS_implicit_ADI
  implicit none
  private
  public :: make_animation_implicit_adi, makeMP4

contains

  subroutine make_animation_implicit_adi(u0, t_end, deltax, deltay, alpha, dt, t1, t2, t3, t4)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha, dt
    real(real64), intent(in) :: t1, t2, t3, t4

    integer :: nmax
    real(real64), allocatable :: u(:,:)
    integer :: every

    write(*,'(/,A)') 'ANIMATION (implicit ADI)'

    nmax = nint(t_end / dt)
    allocate(u(size(u0,1), size(u0,2)))

    every = 1
    call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u, &
                           do_frames=.true., frame_every=every)

    call makeMP4('animation/png/frame_%06d.png', 'animation/movie/implicit_adi.mp4', 30)

    deallocate(u)
  end subroutine make_animation_implicit_adi

  subroutine makeMP4(png_pattern, outmp4, fps)
    character(len=*), intent(in) :: png_pattern
    character(len=*), intent(in) :: outmp4
    integer, intent(in) :: fps

    character(len=2048) :: cmd
    integer :: exitstat

    write(*,'(/,A)') 'Making MP4 with ffmpeg'
    write(*,'(A)')   'Input:  ' // trim(png_pattern)
    write(*,'(A)')   'Output: ' // trim(outmp4)

    write(cmd,'(A,I0,A)') 'ffmpeg -y -framerate ', fps, ' -i "' // trim(png_pattern) // &
                          '" -c:v libx264 -pix_fmt yuv420p "' // trim(outmp4) // '"'

    call execute_command_line(trim(cmd), exitstat=exitstat)

    if (exitstat /= 0) then
      write(*,'(A,I0)') 'ffmpeg failed, exit status = ', exitstat
      write(*,'(A)') 'Check that ffmpeg is installed and the PNG pattern exists.'
    else
      write(*,'(A)') 'MP4 created successfully.'
    end if
  end subroutine makeMP4

end module animation
