! ============================================================================
! 2D Heat Conduction Solver Fortran
!
! Single file, no modules used.
! Every procedure has an explicit interface via CONTAINS.
! ============================================================================

program heat2d
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer :: argc
  character(len=64) :: arg1
  logical :: do_sim, do_conv, do_all, do_anim

  real(real64) :: t_end
  real(real64) :: deltax, deltay, alpha
  real(real64) :: dt_implicit, dt_explicit_given
  real(real64) :: xmin, xmax, ymin, ymax
  integer :: imax, jmax
  real(real64) :: t0, t1, t2, t3, t4

  real(real64), allocatable :: u0(:,:)

  call ensure_project_dirs()

  argc = command_argument_count()
  if (argc >= 1) then
    call get_command_argument(1, arg1)
    arg1 = adjustl(arg1)
  else
    arg1 = ''
  end if

  do_sim  = (argc == 0) .or. (argc >= 1 .and. trim(arg1) == 'sim')
  do_conv = (argc >= 1 .and. trim(arg1) == 'convergence')
  do_all  = (argc >= 1 .and. trim(arg1) == 'all')
  do_anim = (argc >= 1 .and. trim(arg1) == 'animation')

  t_end = 0.5_real64

  deltax = 0.1_real64
  deltay = 0.1_real64
  alpha  = 0.645_real64

  dt_implicit = 0.01_real64
  dt_explicit_given = 0.01_real64

  xmin = 0.0_real64
  xmax = 3.5_real64
  ymin = 0.0_real64
  ymax = 3.5_real64

  imax = int(ceiling((xmax - xmin) / deltax + 1.0_real64))
  jmax = int(ceiling((ymax - ymin) / deltay + 1.0_real64))

  t0 = 0.0_real64
  t1 = 200.0_real64
  t2 = 200.0_real64
  t3 = 0.0_real64
  t4 = 0.0_real64

  allocate(u0(imax, jmax))
  call initializeField(imax, jmax, t0, t1, t2, t3, t4, u0)

  if (do_sim .or. do_all) then
    call sim(u0, t_end, deltax, deltay, alpha, dt_implicit, dt_explicit_given, t1, t2, t3, t4)
  end if

  if (do_conv .or. do_all) then
    call convergence(u0, t_end, deltax, deltay, alpha, t1, t2, t3, t4)
  end if

  if (do_anim) then
    call make_animation_implicit_adi(u0, t_end, deltax, deltay, alpha, dt_implicit, t1, t2, t3, t4)
  end if

  write(*,'(/,A)') 'Done'
  deallocate(u0)

contains

  ! ==========================================================================
  ! Directory helpers
  ! ==========================================================================

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

  ! ==========================================================================
  ! Error norms
  ! ==========================================================================

  real(real64) function error_L2(u, uref) result(val)
    real(real64), intent(in) :: u(:,:), uref(:,:)
    integer :: imax2, jmax2, i, j, n
    real(real64) :: sum, e

    imax2 = size(u,1)
    jmax2 = size(u,2)

    sum = 0.0_real64
    n = 0

    do i = 2, imax2-1
      do j = 2, jmax2-1
        e = u(i,j) - uref(i,j)
        sum = sum + e*e
        n = n + 1
      end do
    end do

    if (n > 0) then
      val = sqrt(sum / real(n, real64))
    else
      val = 0.0_real64
    end if
  end function error_L2

  real(real64) function error_Linf(u, uref) result(val)
    real(real64), intent(in) :: u(:,:), uref(:,:)
    integer :: imax2, jmax2, i, j
    real(real64) :: e

    imax2 = size(u,1)
    jmax2 = size(u,2)

    val = 0.0_real64
    do i = 2, imax2-1
      do j = 2, jmax2-1
        e = abs(u(i,j) - uref(i,j))
        if (e > val) val = e
      end do
    end do
  end function error_Linf

  ! ==========================================================================
  ! High level simulation driver
  ! ==========================================================================

  subroutine sim(u0, t_end, deltax, deltay, alpha, dt_implicit, dt_explicit_given, t1, t2, t3, t4)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, alpha
    real(real64), intent(in) :: dt_implicit, dt_explicit_given
    real(real64), intent(in) :: t1, t2, t3, t4

    integer :: nmax_implicit
    real(real64), allocatable :: u_implicit(:,:), u_explicit(:,:), u_fail(:,:), u_safe(:,:)
    real(real64) :: fx, fy, sum, dt_safe

    write(*,'(/,A)') 'SIMULATING'

    call write_field_xyz('data/initial.dat', u0, deltax, deltay)
    call plotContourMatlabLike('data/initial.dat', 'plot/contour_initial.png', 0.0_real64, 'Initial conditions')

    nmax_implicit = nint(t_end / dt_implicit)
    allocate(u_implicit(size(u0,1), size(u0,2)))

    call FTCS_implicit_ADI(u0, nmax_implicit, deltax, deltay, dt_implicit, alpha, t1, t2, t3, t4, u_implicit)

    call write_field_xyz('data/implicit.dat', u_implicit, deltax, deltay)
    call plotContourMatlabLike('data/implicit.dat', 'plot/contour_implicit.png', t_end, 'Implicit ADI')

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

  ! ==========================================================================
  ! Time convergence study (implicit ADI)
  ! ==========================================================================

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

  ! ==========================================================================
  ! Initialize field
  ! ==========================================================================

  subroutine initializeField(imax, jmax, t0, t1, t2, t3, t4, u)
    integer, intent(in) :: imax, jmax
    real(real64), intent(in) :: t0, t1, t2, t3, t4
    real(real64), intent(out) :: u(imax, jmax)
    integer :: i, j

    u = t0

    do j = 1, jmax
      u(1, j)    = t2
      u(imax, j) = t4
    end do
    do i = 1, imax
      u(i, 1)    = t1
      u(i, jmax) = t3
    end do
  end subroutine initializeField

  ! ==========================================================================
  ! Explicit FTCS solver
  ! ==========================================================================

  subroutine FTCS_Explicit(u0, t_end, deltax, deltay, dt, alpha, t1, t2, t3, t4, u)
    real(real64), intent(in) :: u0(:,:)
    real(real64), intent(in) :: t_end, deltax, deltay, dt, alpha
    real(real64), intent(in) :: t1, t2, t3, t4
    real(real64), intent(out) :: u(size(u0,1), size(u0,2))

    integer :: imax2, jmax2
    integer :: i, j, n, nfull, nsteps
    real(real64) :: dt_last, dt_n, fx, fy
    real(real64), allocatable :: u_new(:,:)

    imax2 = size(u0,1)
    jmax2 = size(u0,2)

    allocate(u_new(imax2, jmax2))
    u = u0

    nfull = int(floor(t_end / dt))
    dt_last = t_end - real(nfull, real64) * dt
    if (dt_last > 0.0_real64) then
      nsteps = nfull + 1
    else
      nsteps = nfull
    end if

    do n = 1, nsteps

      if (n == nsteps .and. dt_last > 0.0_real64) then
        dt_n = dt_last
      else
        dt_n = dt
      end if

      fx = alpha * dt_n / (deltax*deltax)
      fy = alpha * dt_n / (deltay*deltay)

      u_new = u

      do i = 2, imax2-1
        do j = 2, jmax2-1
          u_new(i,j) = (1.0_real64 - 2.0_real64*fx - 2.0_real64*fy) * u(i,j) &
                     + fx * (u(i+1,j) + u(i-1,j)) &
                     + fy * (u(i,j+1) + u(i,j-1))
        end do
      end do

      do j = 1, jmax2
        u_new(1,j)     = t2
        u_new(imax2,j) = t4
      end do
      do i = 1, imax2
        u_new(i,1)     = t1
        u_new(i,jmax2) = t3
      end do

      u = u_new
    end do

    deallocate(u_new)
  end subroutine FTCS_Explicit

  ! ==========================================================================
  ! Write one animation frame
  ! ==========================================================================

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

  ! ==========================================================================
  ! Implicit ADI solver
  ! Optional frame dumping for animation
  ! ==========================================================================

  subroutine FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u, do_frames, frame_every)
    real(real64), intent(in) :: u0(:,:)
    integer, intent(in) :: nmax
    real(real64), intent(in) :: deltax, deltay, dt, alpha
    real(real64), intent(in) :: t1, t2, t3, t4
    real(real64), intent(out) :: u(size(u0,1), size(u0,2))
    logical, intent(in), optional :: do_frames
    integer, intent(in), optional :: frame_every

    integer :: imax2, jmax2
    integer :: i, j, n
    real(real64) :: fx, fy

    real(real64), allocatable :: u_dummy(:,:)

    real(real64), allocatable :: ax(:), bx(:), cx(:), dx(:), solx(:)
    real(real64), allocatable :: ay(:), by(:), cy(:), dy(:), soly(:)

    logical :: frames_on
    integer :: every
    character(len=64) :: scheme

    frames_on = .false.
    if (present(do_frames)) frames_on = do_frames

    every = 1
    if (present(frame_every)) every = max(1, frame_every)

    scheme = 'Implicit ADI animation'

    imax2 = size(u0,1)
    jmax2 = size(u0,2)

    allocate(u_dummy(imax2, jmax2))
    u = u0
    u_dummy = 0.0_real64

    fx = alpha * dt / (deltax*deltax)
    fy = alpha * dt / (deltay*deltay)

    allocate(ax(imax2), bx(imax2), cx(imax2), dx(imax2), solx(imax2))
    ax = 0.0_real64
    bx = 0.0_real64
    cx = 0.0_real64
    dx = 0.0_real64
    solx = 0.0_real64

    do i = 2, imax2-1
      ax(i) = -fx / 2.0_real64
      bx(i) = -fx / 2.0_real64
      dx(i) =  1.0_real64 + fx
    end do
    dx(1) = 1.0_real64
    dx(imax2) = 1.0_real64

    allocate(ay(jmax2), by(jmax2), cy(jmax2), dy(jmax2), soly(jmax2))
    ay = 0.0_real64
    by = 0.0_real64
    cy = 0.0_real64
    dy = 0.0_real64
    soly = 0.0_real64

    do j = 2, jmax2-1
      ay(j) = -fy / 2.0_real64
      by(j) = -fy / 2.0_real64
      dy(j) =  1.0_real64 + fy
    end do
    dy(1) = 1.0_real64
    dy(jmax2) = 1.0_real64

    if (frames_on) then
       call ensure_dir('animation')
       call ensure_dir('animation/dat')
       call ensure_dir('animation/png')

       write(*,'(/,A,I6)') 'Implicit ADI animation: total steps = ', nmax
       write(*,'(A)')     'Writing frames every step'

       write(*,'(A,I6,A,I6)') 'Frame ', 0, ' / ', nmax
       call write_frame(0, u, deltax, deltay, 0.0_real64, scheme)
    end if

    do n = 1, nmax

      do j = 1, jmax2
        u(1,j)     = t2
        u(imax2,j) = t4
      end do
      do i = 1, imax2
        u(i,1)     = t1
        u(i,jmax2) = t3
      end do

      ! X sweep
      do j = 2, jmax2-1

        cx = 0.0_real64
        solx = 0.0_real64

        do i = 2, imax2-1
          cx(i) = (1.0_real64 - fy) * u(i,j) + (fy/2.0_real64) * (u(i,j+1) + u(i,j-1))
        end do

        cx(2)       = cx(2)       + (fx/2.0_real64) * u(1,j)
        cx(imax2-1) = cx(imax2-1) + (fx/2.0_real64) * u(imax2,j)

        solx(1) = u(1,j)
        solx(imax2) = u(imax2,j)

        cx(1) = u(1,j)
        cx(imax2) = u(imax2,j)

        call thomasTriDiagonal(imax2, ax, bx, cx, dx, solx)

        do i = 1, imax2
          u_dummy(i,j) = solx(i)
        end do
      end do

      do j = 1, jmax2
        u_dummy(1,j)     = t2
        u_dummy(imax2,j) = t4
      end do
      do i = 1, imax2
        u_dummy(i,1)     = t1
        u_dummy(i,jmax2) = t3
      end do

      ! Y sweep
      do i = 2, imax2-1

        cy = 0.0_real64
        soly = 0.0_real64

        do j = 2, jmax2-1
          cy(j) = (1.0_real64 - fx) * u_dummy(i,j) + (fx/2.0_real64) * (u_dummy(i+1,j) + u_dummy(i-1,j))
        end do

        cy(2)       = cy(2)       + (fy/2.0_real64) * u_dummy(i,1)
        cy(jmax2-1) = cy(jmax2-1) + (fy/2.0_real64) * u_dummy(i,jmax2)

        soly(1) = u_dummy(i,1)
        soly(jmax2) = u_dummy(i,jmax2)

        cy(1) = soly(1)
        cy(jmax2) = soly(jmax2)

        call thomasTriDiagonal(jmax2, ay, by, cy, dy, soly)

        do j = 1, jmax2
          u(i,j) = soly(j)
        end do
      end do

      do j = 1, jmax2
        u(1,j)     = t2
        u(imax2,j) = t4
      end do
      do i = 1, imax2
        u(i,1)     = t1
        u(i,jmax2) = t3
      end do

      if (frames_on) then
         if (mod(n, every) == 0 .or. n == nmax) then
            write(*,'(A,I6,A,I6)') 'Frame ', n, ' / ', nmax
            call write_frame(n, u, deltax, deltay, real(n, real64)*dt, scheme)
         end if
      end if

    end do

    deallocate(u_dummy, ax, bx, cx, dx, solx, ay, by, cy, dy, soly)
  end subroutine FTCS_implicit_ADI

  ! ==========================================================================
  ! Thomas tridiagonal solver
  ! ==========================================================================

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

  ! ==========================================================================
  ! Write field to file
  ! ==========================================================================

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

  ! ==========================================================================
  ! Plot helpers
  ! ==========================================================================

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

    write(*,'(/,A)') trim(scheme)
    call execute_command_line(trim(cmd))
  end subroutine plotContourMatlabLike

  subroutine plot_time_convergence(infile, outpng)
    character(len=*), intent(in) :: infile, outpng
    character(len=1024) :: cmd

    write(cmd,'(A)') 'gnuplot -e "' // &
                     'datafile=\"' // trim(infile) // '\"; ' // &
                     'outpng=\"'   // trim(outpng) // '\"' // &
                     '" gnuplot_scripts/plot_convergence.gp'

    write(*,'(/,A)') 'Time convergence plot'
    call execute_command_line(trim(cmd))
  end subroutine plot_time_convergence

  ! ==========================================================================
  ! Animation driver (implicit ADI)
  ! ==========================================================================

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
    call FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u, do_frames=.true., frame_every=every)
    call makeMP4('animation/png/frame_%06d.png', 'animation/implicit_adi.mp4', 30)

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

end program heat2d
