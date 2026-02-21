program heat_1d
  use headers, only: dp
  use helpers, only: make_directory
  use driver,  only: sim, convergence_study, plot
  implicit none

  integer  :: imax
  real(dp) :: dx, dt, alpha, t0, tb, Fourier_num
  real(dp) :: t_target

  imax  = 21
  dx    = 0.05_dp
  dt    = 0.01_dp
  alpha = 0.1_dp
  t0    = 100.0_dp
  tb    = 300.0_dp

  Fourier_num = alpha * dt / (dx*dx)
  t_target = 0.4_dp

  call make_directory("data")
  call make_directory("plot")

  call sim(dx, dt, imax, t0, tb, Fourier_num)
  call convergence_study(dx, imax, alpha, t0, tb, t_target)
  call plot(dt, t_target)

  print *, "Done"
end program heat_1d
