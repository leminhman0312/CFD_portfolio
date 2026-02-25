program heat2d
  use kinds,      only: real64
  use io_dirs,    only: ensure_project_dirs
  use field_init, only: initializeField
  use driver,     only: sim, convergence
  use animation,  only: make_animation_implicit_adi
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
  t2 = 0.00_real64
  t3 = 200.00_real64
  t4 = 0.00_real64

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

end program heat2d
