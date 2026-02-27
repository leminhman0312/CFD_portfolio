module driver
  use kinds,       only: real64
  use io_field,    only: write_field_xyz
  use plotting,    only: plotContourMatlabLike, plot_time_convergence
  use solvers_adi, only: FTCS_implicit_ADI
  use norms,       only: error_L2, error_Linf
  implicit none
  private
  public :: sim, convergence

contains

  subroutine sim(u0, t_end, deltax, deltay, alpha, dt_implicit, r_source, t_source)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: dt_implicit
    real(real64), intent(in) :: r_source, t_source

    integer :: nmax
    real(real64), allocatable :: u(:,:)

    write(*,'(/,A)') 'SIMULATING (source case: zero BC + internal hot disk)'

    call write_field_xyz('data/initial.dat', u0, deltax, deltay)
    call plotContourMatlabLike('data/initial.dat', 'plot/contour_initial.png', 0.0_real64, 'Initial conditions')

    nmax = nint(t_end / dt_implicit)
    allocate(u(size(u0,1), size(u0,2)))

    call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt_implicit, alpha, &
                           0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u, &
                           use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)

    call write_field_xyz('data/implicit_source.dat', u, deltax, deltay)
    call plotContourMatlabLike('data/implicit_source.dat', 'plot/contour_implicit_source.png', t_end, 'Implicit ADI source')

    deallocate(u)
  end subroutine sim


  subroutine convergence(u0, t_end, deltax, deltay, alpha, r_source, t_source)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: r_source, t_source

    real(real64) :: dt_start, dt_ratio, dt_min
    real(real64), allocatable :: dt_list(:)
    real(real64) :: dt, dt_ref
    integer :: n_dt, k, nmax, nmax_ref, unit
    real(real64), allocatable :: u_ref(:,:), u(:,:)
    real(real64) :: L2, Linf

    write(*,'(/,A)') 'TIME CONVERGENCE STUDY (implicit ADI, source case)'
    write(*,'(A)') 'dt        L2 error      Linf error'
    write(*,'(A)') '----------------------------------'

    dt_start = 0.04_real64
    dt_ratio = 0.5_real64
    dt_min   = 1.0e-6_real64

    n_dt = 0
    dt = dt_start
    do while (dt >= dt_min)
      n_dt = n_dt + 1
      dt = dt * dt_ratio
    end do

    allocate(dt_list(n_dt))
    dt = dt_start
    do k = 1, n_dt
      dt_list(k) = dt
      dt = dt * dt_ratio
    end do

    dt_ref = dt_list(n_dt)
    nmax_ref = nint(t_end / dt_ref)

    allocate(u_ref(size(u0,1), size(u0,2)))
    call FTCS_implicit_ADI(u0, nmax_ref, deltax, deltay, dt_ref, alpha, &
                           0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u_ref, &
                           use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)

    open(newunit=unit, file='data/time_convergence.dat', status='replace', action='write')

    allocate(u(size(u0,1), size(u0,2)))

    do k = 1, n_dt
      dt = dt_list(k)
      nmax = nint(t_end / dt)

      call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, &
                             0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u, &
                             use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)

      L2   = error_L2(u, u_ref)
      Linf = error_Linf(u, u_ref)

      write(*,'(F8.5,2X,ES12.6,2X,ES12.6)') dt, L2, Linf
      write(unit,'(F12.10,1X,ES16.10,1X,ES16.10)') dt, L2, Linf
    end do

    close(unit)

    write(*,'(/,A)') 'Wrote data/time_convergence.dat'
    call plot_time_convergence('data/time_convergence.dat', 'plot/time_convergence.png')

    deallocate(dt_list, u_ref, u)
  end subroutine convergence

end module driver
