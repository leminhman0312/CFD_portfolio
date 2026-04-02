program heat2d
  use kinds,      only: real64
  use io_dirs,    only: ensure_project_dirs
  use field_init, only: initializeField
  use driver,     only: sim, convergence
  implicit none

  integer :: argc
  character(len=64) :: arg1
  character(len=16) :: cmd

  logical :: do_sim, do_conv, do_anim, do_all

  real(real64) :: t_end
  real(real64) :: deltax, deltay, alpha
  real(real64) :: dt_implicit, dt_explicit_given
  real(real64) :: xmin, xmax, ymin, ymax
  integer :: imax, jmax
  real(real64) :: t0, t1, t2, t3, t4

  real(real64), allocatable :: u0(:,:)

  call ensure_project_dirs()

  argc = command_argument_count()
  arg1 = ''
  if (argc >= 1) then
    call get_command_argument(1, arg1)
    arg1 = adjustl(arg1)
  end if

  if (argc == 0) then
    cmd = 'sim'
  else
    cmd = trim(arg1)
  end if

  do_sim  = (cmd == 'sim')        .or. (cmd == 'all')
  do_conv = (cmd == 'convergence').or. (cmd == 'all')
  do_anim = (cmd == 'animation') .or. (cmd == 'all')
  do_all  = (cmd == 'all')

  t_end = 2.00_real64

  deltax = 0.1_real64
  deltay = 0.1_real64
  alpha  = 0.645_real64

  dt_implicit       = 0.01_real64
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

  if (do_sim .or. do_anim) then
    call sim(u0, t_end, deltax, deltay, alpha, &
             dt_implicit, dt_explicit_given, &
             t1, t2, t3, t4, anim_on=do_anim)
  end if

  if (do_conv) then
    call convergence(u0, t_end, deltax, deltay, alpha, t1, t2, t3, t4)
  end if

  write(*,'(/,A)') 'Done'
  deallocate(u0)

end program heat2d
