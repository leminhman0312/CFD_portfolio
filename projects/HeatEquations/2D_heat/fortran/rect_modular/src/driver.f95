module driver
  use kinds,            only: real64
  use io_field,         only: write_field_xyz
  use plotting,         only: plotContourMatlabLike, plot_time_convergence
  use solvers_explicit, only: FTCS_Explicit
  use solvers_adi,      only: FTCS_implicit_ADI
  use norms,            only: error_L2, error_Linf
  use animation,  only: makeMP4, clear_animation_frames
  implicit none
  private
  public :: sim, convergence

contains

    subroutine sim(u0, t_end, deltax, deltay, alpha, dt_implicit, dt_explicit_given, t1, t2, t3, t4, anim_on)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: dt_implicit, dt_explicit_given
    real(real64), intent(in) :: t1, t2, t3, t4
    logical, intent(in), optional :: anim_on

    integer :: nmax_implicit
    real(real64), allocatable :: u_implicit(:,:), u_explicit(:,:), u_fail(:,:), u_safe(:,:)
    real(real64) :: fx, fy, sum, dt_safe
    logical :: do_anim

    do_anim = .false.
    if (present(anim_on)) do_anim = anim_on

    write(*,'(/,A)') 'SIMULATING'

    call write_field_xyz('data/initial.dat', u0, deltax, deltay)
    call plotContourMatlabLike('data/initial.dat', 'plot/contour_initial.png', 0.0_real64, 'Initial conditions')

    ! ---------------------------
    ! Implicit ADI solution
    ! ---------------------------
    nmax_implicit = nint(t_end / dt_implicit)
    allocate(u_implicit(size(u0,1), size(u0,2)))

    write(*,*) 'Running FTCS implicit'
    if (do_anim) then
      call clear_animation_frames()
    end if

    call FTCS_implicit_ADI(u0, nmax_implicit, deltax, deltay, dt_implicit, alpha, t1, t2, t3, t4, u_implicit, &
                           do_frames=do_anim, frame_every=1)

    call write_field_xyz('data/implicit.dat', u_implicit, deltax, deltay)
    call plotContourMatlabLike('data/implicit.dat', 'plot/contour_implicit.png', t_end, 'Implicit ADI')

    if (do_anim) then
      call makeMP4('animation/png/frame_%06d.png', 'animation/movie/implicit_adi.mp4', 30)
    end if

    ! ---------------------------
    ! Explicit FTCS stability logic
    ! ---------------------------

    write(*,*) 'Running FTCS Explicit'
    fx = alpha * dt_explicit_given / (deltax*deltax)
    fy = alpha * dt_explicit_given / (deltay*deltay)
    sum = fx + fy

    if (sum <= 0.5_real64) then

      write(*,*) 'FTCS Explicit: fx + fy = ', sum, ' <= 0.5'
      allocate(u_explicit(size(u0,1), size(u0,2)))

      if (do_anim) then
        call clear_animation_frames()
      end if

      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_explicit_given, alpha, t1, t2, t3, t4, u_explicit, &
                         anim_on=do_anim, frame_every=1)

      call write_field_xyz('data/explicit.dat', u_explicit, deltax, deltay)
      call plotContourMatlabLike('data/explicit.dat', 'plot/contour_explicit.png', t_end, 'Explicit FTCS')

      if (do_anim) then
        call makeMP4('animation/png/frame_%06d.png', 'animation/movie/explicit.mp4', 30)
      end if

      deallocate(u_explicit)

   else

       write(*,*) 'FTCS Explicit: fx + fy = ', sum, ' > 0.5'
       write(*,*) 'Calculating safe dt'
      dt_safe = 0.5_real64 / (alpha * (1.0_real64/(deltax*deltax) + 1.0_real64/(deltay*deltay)))

      allocate(u_fail(size(u0,1), size(u0,2)))
      allocate(u_safe(size(u0,1), size(u0,2)))

      ! Run the unstable dt case (animate if requested)

      write(*,*) 'Running unstable dt'
      if (do_anim) then
        call clear_animation_frames()
      end if

      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_explicit_given, alpha, t1, t2, t3, t4, u_fail, &
                         anim_on=do_anim, frame_every=1)

      call write_field_xyz('data/explicit_failed.dat', u_fail, deltax, deltay)
      call plotContourMatlabLike('data/explicit_failed.dat', 'plot/contour_explicit_failed.png', t_end, 'Explicit unstable dt')

      if (do_anim) then
        call makeMP4('animation/png/frame_%06d.png', 'animation/movie/explicit_failed.mp4', 30)
      end if

      ! Run the safe dt case (animate if requested)

      write(*,*) 'Running safe dt'
      if (do_anim) then
        call clear_animation_frames()
      end if

      call FTCS_Explicit(u0, t_end, deltax, deltay, dt_safe, alpha, t1, t2, t3, t4, u_safe, &
                         anim_on=do_anim, frame_every=1)

      call write_field_xyz('data/explicit_safe.dat', u_safe, deltax, deltay)
      call plotContourMatlabLike('data/explicit_safe.dat', 'plot/contour_explicit_safe.png', t_end, 'Explicit safe dt')

      if (do_anim) then
        call makeMP4('animation/png/frame_%06d.png', 'animation/movie/explicit_safe.mp4', 30)
      end if

      deallocate(u_fail, u_safe)
    end if

    deallocate(u_implicit)
  end subroutine sim



  subroutine convergence(u0, t_end, deltax, deltay, alpha, t1, t2, t3, t4)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: t1, t2, t3, t4

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

    call FTCS_implicit_ADI(u0, nmax_ref, deltax, deltay, dt_ref, alpha, t1, t2, t3, t4, u_ref)

    open(newunit=unit, file='data/time_convergence.dat', status='replace', action='write')

    allocate(u(size(u0,1), size(u0,2)))

    do k = 1, n_dt
      dt = dt_list(k)
      nmax = nint(t_end / dt)

      call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u)

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
