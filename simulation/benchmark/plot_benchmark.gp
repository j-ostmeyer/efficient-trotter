set terminal postscript landscape enhanced color lw 1 font 16
set encoding iso_8859_1
set grid lt 0 lw .3
set samples 10000
set fit results
set pointsize 1
set border lw 0.5

set linetype 1 pt 4
set linetype 2 pt 6
set linetype 3 pt 8
set linetype 4 pt 10
set linetype 5 pt 1
set linetype 6 pt 2
set linetype 7 pt 2

!rm -f *.eps *.tex

schemes = "Leap-Frog Omelyan2 Forest-Ruth FR-Type small-B Non-Unitary1 Suzuki4 Algorithm-30 Optimal-4th-order Non-Unitary2 Omelyan-ST-4 Non-Unitary3 Non-Unit-Blanes Uniform-Non-Unit Blanes-4 Symplectic-6 FR-Squared Blanes-6 FR-Suzuki6 Suzuki6 FR-Cubed FR-Suzuki8 Blanes6-Suzuki8 Suzuki8 Taylor"
array cyc[25] = [1, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 9, 10, 15, 25, 27, 45, 50, 125, 3]
array algs[19] = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 18, 20, 23, 24, 25]

L = 6

set logscale xy
set key bottom left
set xlabel "cycles/dt"
set ylabel "error"

set output "2-stage_fix-t.eps"
set title "2 Stages, fixed t"
plot [2:][:1] for [i=1:19] "2-stage_fix-t.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])

set output "2L-stage_fix-t.eps"
set title "2L Stages, fixed t"
plot [2:][:1] for [i=1:19] "2L-stage_fix-t.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])

set output "3-stage_fix-t.eps"
set title "3 Stages, fixed t"
plot [2:][:1] for [i=1:19] "3-stage_fix-t.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])

set output "3L-stage_fix-t.eps"
set title "3L Stages, fixed t"
plot [2:][:1] for [i=1:19] "3L-stage_fix-t.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])

#set output "2-stage_fix-t_small-x.eps"
#set title "2 Stages, fixed t, small X"
#plot [2:][:1] for [i=1:19] "2-stage_fix-t_small-x.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])
#
#set output "2L-stage_fix-t_small-x.eps"
#set title "2L Stages, fixed t, small X"
#plot [2:][:1] for [i=1:19] "2L-stage_fix-t_small-x.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])
#
#set output "3-stage_fix-t_small-x.eps"
#set title "3 Stages, fixed t, small XY"
#plot [2:][:1] for [i=1:19] "3-stage_fix-t_small-x.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])
#
#set output "3L-stage_fix-t_small-x.eps"
#set title "3L Stages, fixed t, small XY"
#plot [2:][:1] for [i=1:19] "3L-stage_fix-t_small-x.txt" u 1:(column(algs[i]+1)/2**L) w l dt i t word(schemes, algs[i])

unset logscale x
set key bottom right
set xlabel "t"
set ylabel "error/t"

set output "2-stage_fix-dt.eps"
set title "2 Stages, fixed dt"
plot [:10] for [i=1:19] "2-stage_fix-dt.txt" u ($1*cyc[algs[i]]):(column(algs[i]+1)/2**L/($1*cyc[algs[i]])) w l dt i t word(schemes, algs[i])

set output "2L-stage_fix-dt.eps"
set title "2L Stages, fixed dt"
plot [:10] for [i=1:19] "2L-stage_fix-dt.txt" u ($1*cyc[algs[i]]):(column(algs[i]+1)/2**L/($1*cyc[algs[i]])) w l dt i t word(schemes, algs[i])

set output "3-stage_fix-dt.eps"
set title "3 Stages, fixed dt"
plot [:10] for [i=1:19] "3-stage_fix-dt.txt" u ($1*cyc[algs[i]]):(column(algs[i]+1)/2**L/($1*cyc[algs[i]])) w l dt i t word(schemes, algs[i])

set output "3L-stage_fix-dt.eps"
set title "3L Stages, fixed dt"
plot [:10] for [i=1:19] "3L-stage_fix-dt.txt" u ($1*cyc[algs[i]]):(column(algs[i]+1)/2**L/($1*cyc[algs[i]])) w l dt i t word(schemes, algs[i])
