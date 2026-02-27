module driver
  use kinds,            only: real64
  use io_field,         only: write_field_xyz
  use plotting,         only: plotContourMatlabLike, plot_time_convergence
  use solvers_explicit, only: FTCS_Explicit
  use solvers_adi,      only: FTCS_implicit_ADI
  use norms,            only: error_L2, error_Linf
  implicit none
  private
  public :: sim, convergence

contains

  subroutine sim(u0, t_end, deltax, deltay, alpha, dt_implicit, dt_explicit_given, t1, t2, t3, t4, &
                 bc_case, r_source, t_source)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: dt_implicit, dt_explicit_given
    real(real64), intent(in) :: t1, t2, t3, t4
    character(len=*), intent(in) :: bc_case
    real(real64), intent(in) :: r_source, t_source

    integer :: nmax_implicit
    real(real64), allocatable :: u_implicit(:,:), u_explicit(:,:), u_fail(:,:), u_safe(:,:)
    real(real64) :: fx, fy, sum, dt_safe

    write(*,'(/,A)') 'SIMULATING'

    call write_field_xyz('data/initial.dat', u0, deltax, deltay)
    call plotContourMatlabLike('data/initial.dat', 'plot/contour_initial.png', 0.0_real64, 'Initial conditions')

    nmax_implicit = nint(t_end / dt_implicit)
    allocate(u_implicit(size(u0,1), size(u0,2)))

    if (trim(bc_case) == 'source') then
      call FTCS_implicit_ADI(u0, nmax_implicit, deltax, deltay, dt_implicit, alpha, &
                             0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u_implicit, &
                             use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)
      call write_field_xyz('data/implicit_source.dat', u_implicit, deltax, deltay)
      call plotContourMatlabLike('data/implicit_source.dat', 'plot/contour_implicit_source.png', t_end, 'Implicit ADI source')
    else
      call FTCS_implicit_ADI(u0, nmax_implicit, deltax, deltay, dt_implicit, alpha, t1, t2, t3, t4, u_implicit)
      call write_field_xyz('data/implicit.dat', u_implicit, deltax, deltay)
      call plotContourMatlabLike('data/implicit.dat', 'plot/contour_implicit.png', t_end, 'Implicit ADI')
    end if

    if (trim(bc_case) == 'source') then
      write(*,'(A)') 'Skipping explicit FTCS for source case.'
      deallocate(u_implicit)
      return
    end if

    fx = alpha * dt_explicit_given / (deltax*deltax)
    fy = alpha * dt_explicit_given / (deltay*deltay)
    sum = fx + fy

    if (sum <= 0.5_real64) then
      allocate(u_explicit(size(u0,1), size(u0,2)))
      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_explicit_given, alpha, t1, t2, t3, t4, u_explicit)
      call write_field_xyz('data/explicit.dat', u_explicit, deltax, deltay)
      call plotContourMatlabLike('data/explicit.dat', 'plot/contour_explicit.png', t_end, 'Explicit FTCS')
      deallocate(u_explicit)
    else
      dt_safe = 0.5_real64 / (alpha * (1.0_real64/(deltax*deltax) + 1.0_real64/(deltay*deltay)))

      allocate(u_fail(size(u0,1), size(u0,2)))
      allocate(u_safe(size(u0,1), size(u0,2)))

      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_explicit_given, alpha, t1, t2, t3, t4, u_fail)
      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_safe,           alpha, t1, t2, t3, t4, u_safe)

      call write_field_xyz('data/explicit_failed.dat', u_fail, deltax, deltay)
      call write_field_xyz('data/explicit_safe.dat',   u_safe, deltax, deltay)

      call plotContourMatlabLike('data/explicit_failed.dat', 'plot/contour_explicit_failed.png', t_end, 'Explicit unstable dt')
      call plotContourMatlabLike('data/explicit_safe.dat',   'plot/contour_explicit_safe.png',   t_end, 'Explicit safe dt')

      deallocate(u_fail, u_safe)
    end if

    deallocate(u_implicit)
  end subroutine sim


  subroutine convergence(u0, t_end, deltax, deltay, alpha, t1, t2, t3, t4, &
                         bc_case, r_source, t_source)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: t1, t2, t3, t4
    character(len=*), intent(in) :: bc_case
    real(real64), intent(in) :: r_source, t_source

    real(real64) :: dt_start, dt_ratio, dt_min
    real(real64), allocatable :: dt_list(:)
    real(real64) :: dt, dt_ref
    integer :: n_dt, k, nmax, nmax_ref, unit
    real(real64), allocatable :: u_ref(:,:), u(:,:)
    real(real64) :: L2, Linf

    write(*,'(/,A)') 'TIME CONVERGENCE STUDY (implicit ADI)'
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

    if (trim(bc_case) == 'source') then
      call FTCS_implicit_ADI(u0, nmax_ref, deltax, deltay, dt_ref, alpha, &
                             0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u_ref, &
                             use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)
    else
      call FTCS_implicit_ADI(u0, nmax_ref, deltax, deltay, dt_ref, alpha, t1, t2, t3, t4, u_ref)
    end if

    open(newunit=unit, file='data/time_convergence.dat', status='replace', action='write')

    allocate(u(size(u0,1), size(u0,2)))

    do k = 1, n_dt
      dt = dt_list(k)
      nmax = nint(t_end / dt)

      if (trim(bc_case) == 'source') then
        call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, &
                               0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, u, &
                               use_zero_bc=.true., use_source=.true., r_source=r_source, t_source=t_source)
      else
        call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u)
      end if

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
