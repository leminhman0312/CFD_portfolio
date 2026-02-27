# Time convergence plot for 2D heat solver

set terminal pngcairo size 900,700 enhanced font "Helvetica,12"
set output outpng

set title "Time Convergence Study (Implicit ADI)"
set xlabel "Time step dt"
set ylabel "Error norm"

set logscale x
set logscale y

# Reverse x-axis so dt decreases left -> right
set xrange [*:*] reverse

set grid
set key top right

plot \
  datafile using 1:2 with linespoints lw 2 pt 7 title "L2 error", \
  datafile using 1:3 with linespoints lw 2 pt 9 title "Linf error"
