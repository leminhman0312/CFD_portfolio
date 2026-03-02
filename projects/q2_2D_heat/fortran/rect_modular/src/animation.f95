module animation
  use kinds,     only: real64
  use io_dirs,   only: ensure_animation_dirs
  use io_field,  only: write_field_xyz
  use plotting,  only: plotContourMatlabLike
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  private
  public :: makeMP4, clear_animation_frames, write_frame, animation_step

contains

  subroutine clear_animation_frames()
    integer :: exitstat
    call ensure_animation_dirs()
    call execute_command_line('rm -f -- animation/dat/*.dat animation/png/*.png', exitstat=exitstat)
  end subroutine clear_animation_frames


  subroutine animation_step(step, nmax, scheme, anim_on, every)
    integer, intent(in) :: step, nmax
    character(len=*), intent(in) :: scheme
    logical, intent(in) :: anim_on
    integer, intent(in) :: every

    if (.not. anim_on) return
    if (mod(step, every) /= 0 .and. step /= nmax) return

    write(output_unit,'(A,A,I6,A,I6)', advance='no') &
         char(13), trim(scheme)//': Frame ', step, ' / ', nmax
    flush(output_unit)
  end subroutine animation_step


  subroutine write_frame(step, u, deltax, deltay, time_hr, scheme)
    integer, intent(in) :: step
    real(real64), intent(in) :: u(:,:)
    real(real64), intent(in) :: deltax, deltay, time_hr
    character(len=*), intent(in) :: scheme

    character(len=256) :: datfile, pngfile

    write(datfile,'(A,I6.6,A)') 'animation/dat/frame_', step, '.dat'
    write(pngfile,'(A,I6.6,A)') 'animation/png/frame_', step, '.png'

    call write_field_xyz(trim(datfile), u, deltax, deltay)
    call plotContourMatlabLike(trim(datfile), trim(pngfile), time_hr, scheme)
  end subroutine write_frame


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
