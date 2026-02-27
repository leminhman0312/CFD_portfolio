module solvers_adi
  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds,     only: real64
  use tridiag,   only: thomasTriDiagonal
  use io_dirs,   only: ensure_animation_dirs
  use io_field,  only: write_field_xyz
  use plotting,  only: plotContourMatlabLike
  implicit none
  private
  public :: FTCS_implicit_ADI

contains

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


  subroutine apply_bc_dirichlet(u, t1, t2, t3, t4)
    real(real64), intent(inout) :: u(:,:)
    real(real64), intent(in) :: t1, t2, t3, t4
    integer :: im, jm, i, j

    im = size(u,1)
    jm = size(u,2)

    do j = 1, jm
      u(1, j)  = t2
      u(im, j) = t4
    end do
    do i = 1, im
      u(i, 1)  = t1
      u(i, jm) = t3
    end do
  end subroutine apply_bc_dirichlet


  subroutine apply_bc_zero(u)
    real(real64), intent(inout) :: u(:,:)
    integer :: im, jm, i, j

    im = size(u,1)
    jm = size(u,2)

    do j = 1, jm
      u(1, j)  = 0.0_real64
      u(im, j) = 0.0_real64
    end do
    do i = 1, im
      u(i, 1)  = 0.0_real64
      u(i, jm) = 0.0_real64
    end do
  end subroutine apply_bc_zero


  subroutine clamp_circle(u, deltax, deltay, r, tval)
    real(real64), intent(inout) :: u(:,:)
    real(real64), intent(in) :: deltax, deltay, r, tval

    integer :: im, jm, i, j
    real(real64) :: x0, y0, x, y, dx, dy, rr2

    im = size(u,1)
    jm = size(u,2)

    x0 = 0.5_real64 * real(im - 1, real64) * deltax
    y0 = 0.5_real64 * real(jm - 1, real64) * deltay
    rr2 = r * r

    do j = 1, jm
      y  = real(j - 1, real64) * deltay
      dy = y - y0
      do i = 1, im
        x  = real(i - 1, real64) * deltax
        dx = x - x0
        if (dx*dx + dy*dy <= rr2) then
          u(i,j) = tval
        end if
      end do
    end do
  end subroutine clamp_circle


  subroutine enforce_bc_and_source(u, deltax, deltay, zbc, src, t1, t2, t3, t4, rs, ts)
    real(real64), intent(inout) :: u(:,:)
    real(real64), intent(in) :: deltax, deltay
    logical, intent(in) :: zbc, src
    real(real64), intent(in) :: t1, t2, t3, t4
    real(real64), intent(in) :: rs, ts

    if (zbc) then
      call apply_bc_zero(u)
    else
      call apply_bc_dirichlet(u, t1, t2, t3, t4)
    end if

    if (src) then
      call clamp_circle(u, deltax, deltay, rs, ts)
    end if
  end subroutine enforce_bc_and_source


  subroutine FTCS_implicit_ADI(u0, nmax, deltax, deltay, dt, alpha, t1, t2, t3, t4, u, &
                               do_frames, frame_every, use_zero_bc, use_source, r_source, t_source)
    real(real64), intent(in) :: u0(:,:)
    integer, intent(in) :: nmax
    real(real64), intent(in) :: deltax, deltay, dt, alpha
    real(real64), intent(in) :: t1, t2, t3, t4
    real(real64), intent(out) :: u(size(u0,1), size(u0,2))
    logical, intent(in), optional :: do_frames
    integer, intent(in), optional :: frame_every
    logical, intent(in), optional :: use_zero_bc
    logical, intent(in), optional :: use_source
    real(real64), intent(in), optional :: r_source, t_source

    integer :: imax2, jmax2
    integer :: i, j, n
    real(real64) :: fx, fy

    real(real64), allocatable :: u_dummy(:,:)
    real(real64), allocatable :: ax(:), bx(:), cx(:), dx(:), solx(:)
    real(real64), allocatable :: ay(:), by(:), cy(:), dy(:), soly(:)

    logical :: frames_on
    integer :: every
    character(len=64) :: scheme

    logical :: zbc, src
    real(real64) :: rs, ts

    frames_on = .false.
    if (present(do_frames)) frames_on = do_frames

    every = 1
    if (present(frame_every)) every = max(1, frame_every)

    zbc = .false.
    src = .false.
    rs  = 0.0_real64
    ts  = 0.0_real64

    if (present(use_zero_bc)) zbc = use_zero_bc
    if (present(use_source))  src = use_source
    if (present(r_source))    rs  = r_source
    if (present(t_source))    ts  = t_source

    scheme = 'Implicit ADI'

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

    call enforce_bc_and_source(u, deltax, deltay, zbc, src, t1, t2, t3, t4, rs, ts)

    if (frames_on) then
      call ensure_animation_dirs()
      write(*,'(/,A,I0)') 'Implicit ADI animation: total steps = ', nmax
      write(*,'(A,I0)') 'Writing frames every ', every
      write(*,'(A,I0,A,I0)', advance='no') 'Frame ', 0, ' / ', nmax
      flush(output_unit)
      call write_frame(0, u, deltax, deltay, 0.0_real64, scheme)
    end if

    do n = 1, nmax

      call enforce_bc_and_source(u, deltax, deltay, zbc, src, t1, t2, t3, t4, rs, ts)

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

      call enforce_bc_and_source(u_dummy, deltax, deltay, zbc, src, t1, t2, t3, t4, rs, ts)

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

      call enforce_bc_and_source(u, deltax, deltay, zbc, src, t1, t2, t3, t4, rs, ts)

      if (frames_on) then
        if (mod(n, every) == 0 .or. n == nmax) then
          write(*,'(A)', advance='no') char(13)
          write(*,'(A,I0,A,I0)', advance='no') 'Frame ', n, ' / ', nmax
          flush(output_unit)
          call write_frame(n, u, deltax, deltay, real(n, real64)*dt, scheme)
        end if
      end if

    end do

    if (frames_on) then
      write(*,*)
    end if

    deallocate(u_dummy, ax, bx, cx, dx, solx, ay, by, cy, dy, soly)
  end subroutine FTCS_implicit_ADI

end module solvers_adi
