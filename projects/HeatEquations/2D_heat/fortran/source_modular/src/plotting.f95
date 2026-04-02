module plotting
  use kinds, only: real64
  implicit none
  private
  public :: plotContourMatlabLike, plot_time_convergence

contains

  subroutine plotContourMatlabLike(datafile, outpng, time_hr, scheme)
    character(len=*), intent(in) :: datafile, outpng, scheme
    real(real64), intent(in) :: time_hr

    character(len=1024) :: cmd
    character(len=64) :: tbuf

    write(tbuf,'(F0.3)') time_hr

    write(cmd,'(A)') 'gnuplot -e "' // &
                     'datafile=\"' // trim(datafile) // '\"; ' // &
                     'outpng=\"'   // trim(outpng)   // '\"; ' // &
                     'tlabel='     // trim(adjustl(tbuf)) // '; ' // &
                     'scheme=\"'   // trim(scheme) // '\"' // &
                     '" gnuplot_scripts/plot_contour_2d_matlab_like.gp'

    call execute_command_line(trim(cmd))
  end subroutine plotContourMatlabLike

  subroutine plot_time_convergence(infile, outpng)
    character(len=*), intent(in) :: infile, outpng
    character(len=1024) :: cmd

    write(cmd,'(A)') 'gnuplot -e "' // &
                     'datafile=\"' // trim(infile) // '\"; ' // &
                     'outpng=\"'   // trim(outpng) // '\"' // &
                     '" gnuplot_scripts/plot_convergence.gp'

    call execute_command_line(trim(cmd))
  end subroutine plot_time_convergence

end module plotting
