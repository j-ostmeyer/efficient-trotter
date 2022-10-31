set terminal epslatex color lw 4
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

!rm -f *.eps *.tex

schemes = "Verlet~$(2,1)$~\\eqref{eq:1-leap-frog} Omelyan~$(2,2)$~\\eqref{eq:2-omelyan2} Forest-Ruth~$(4,3)$~\\eqref{eq:3-forest-ruth} FR-Type~$(4,4)$~\\eqref{eq:4-fr-type} small~$A$~$(4,4)$~\\eqref{eq:4-small-B} Non-Unitary~$(4,4)$~\\eqref{eq:4-non-unitary1} Suzuki~$(4,5)$~\\eqref{eq:5-suzuki4} Algorithm-30 Opt.~4th~order~$(4,5)$~\\eqref{eq:5-opt-4th-ord} Non-Unitary~$(4,5)$~\\eqref{eq:5-non-unitary2} Omelyan-ST-4 Non-Unitary3 Non-Unitary-Blanes Unif.~Non-Unitary~$(4,5)$~\\eqref{eq:5-non-unitary-const} Blanes\\&Moan~$(4,6)$~\\eqref{eq:6-blanes4} Symplectic-6 FR-Squared Blanes\\&Moan~$(6,10)$~\\eqref{eq:10-blanes6} FR-Suzuki6 Suzuki~$(6,25)$~(sec.~\\ref{sec:suzuki6}) FR-Cubed FR-Suzuki8 BM6+S~$(8,50)$~(sec.~\\ref{sec:very_high_orders}) Suzuki~$(8,125)$~(sec.~\\ref{sec:very_high_orders}) Taylor~(sec.~\\ref{sec:taylor_expansion})"
array cyc[25] = [1, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 9, 10, 15, 25, 27, 45, 50, 125, 3]

array algs4[8] = [3, 4, 6, 7, 9, 10, 14, 15]
array algsRest[7] = [1, 2, 18, 20, 23, 24, 25]
array algsOpt[7] = [1, 6, 7, 15, 18, 23, 25]

L = 6

set format y "$10^{%T}$"
set logscale xy
set key bottom left
set key width -25
set xlabel "$q/h$"
set ylabel "error"

set output "2-stage_fix-t_ord4.tex"
plot [2:50][1e-6:2] for [i=1:8] "2-stage_fix-t.txt" u 1:(column(algs4[i]+1)/2**(L/2)) w l dt i t word(schemes, algs4[i])

set output "2L-stage_fix-t_ord4.tex"
plot [2:50][1e-6:2] for [i=1:8] "2L-stage_fix-t.txt" u 1:(column(algs4[i]+1)/2**(L/2)) w l dt i t word(schemes, algs4[i])

set output "3-stage_fix-t_ord4.tex"
plot [2:50][1e-6:2] for [i=1:8] "3-stage_fix-t.txt" u 1:(column(algs4[i]+1)/2**(L/2)) w l dt i t word(schemes, algs4[i])

set output "3L-stage_fix-t_ord4.tex"
plot [2:50][1e-6:2] for [i=1:8] "3L-stage_fix-t.txt" u 1:(column(algs4[i]+1)/2**(L/2)) w l dt i t word(schemes, algs4[i])

unset logscale x
set key top right
set xlabel "$t$"
set ylabel "$\\mathrm{error}/t$"

set output "2-stage_fix-dt_ord4.tex"
plot [:10][2e-8:1.5e-5] for [i=1:8] "2-stage_fix-dt.txt" u ($1*cyc[algs4[i]]):(column(algs4[i]+1)/2**(L/2)/($1*cyc[algs4[i]])) w l dt i t word(schemes, algs4[i])

#set output "2L-stage_fix-dt_ord4.tex"
#plot [:10] for [i=1:7] "2L-stage_fix-dt.txt" u ($1*cyc[algs4[i]]):(column(algs4[i]+1)/2**(L/2)/($1*cyc[algs4[i]])) w l dt i t word(schemes, algs4[i])
#
#set output "3-stage_fix-dt_ord4.tex"
#plot [:10] for [i=1:7] "3-stage_fix-dt.txt" u ($1*cyc[algs4[i]]):(column(algs4[i]+1)/2**(L/2)/($1*cyc[algs4[i]])) w l dt i t word(schemes, algs4[i])

set key bottom right
set key height 2
set output "3L-stage_fix-dt_ord4.tex"
plot [:10][6e-8:5e-5] for [i=1:8] "3L-stage_fix-dt.txt" u ($1*cyc[algs4[i]]):(column(algs4[i]+1)/2**(L/2)/($1*cyc[algs4[i]])) w l dt i t word(schemes, algs4[i])

set logscale xy
set key width -18
set key bottom left
set key height 3
set xlabel "$q/h$"
set ylabel "error"

set output "2-stage_fix-t_ordRest.tex"
plot [2:][:3] for [i=1:7] "2-stage_fix-t.txt" u 1:(column(algsRest[i]+1)/2**(L/2)) w l dt i t word(schemes, algsRest[i])

set output "2L-stage_fix-t_ordRest.tex"
plot [2:][:3] for [i=1:7] "2L-stage_fix-t.txt" u 1:(column(algsRest[i]+1)/2**(L/2)) w l dt i t word(schemes, algsRest[i])

set output "3-stage_fix-t_ordRest.tex"
plot [2:][:3] for [i=1:7] "3-stage_fix-t.txt" u 1:(column(algsRest[i]+1)/2**(L/2)) w l dt i t word(schemes, algsRest[i])

set output "3L-stage_fix-t_ordRest.tex"
plot [2:][:3] for [i=1:7] "3L-stage_fix-t.txt" u 1:(column(algsRest[i]+1)/2**(L/2)) w l dt i t word(schemes, algsRest[i])

unset logscale x
set key bottom right
set xlabel "$t$"
set ylabel "$\\mathrm{error}/t$"

set output "2-stage_fix-dt_ordRest.tex"
plot [:10] for [i=1:7] "2-stage_fix-dt.txt" u ($1*cyc[algsRest[i]]):(column(algsRest[i]+1)/2**(L/2)/($1*cyc[algsRest[i]])) w l dt i t word(schemes, algsRest[i])

#set output "2L-stage_fix-dt_ordRest.tex"
#plot [:10] for [i=1:7] "2L-stage_fix-dt.txt" u ($1*cyc[algsRest[i]]):(column(algsRest[i]+1)/2**(L/2)/($1*cyc[algsRest[i]])) w l dt i t word(schemes, algsRest[i])
#
#set output "3-stage_fix-dt_ordRest.tex"
#plot [:10] for [i=1:7] "3-stage_fix-dt.txt" u ($1*cyc[algsRest[i]]):(column(algsRest[i]+1)/2**(L/2)/($1*cyc[algsRest[i]])) w l dt i t word(schemes, algsRest[i])

set output "3L-stage_fix-dt_ordRest.tex"
plot [:10] for [i=1:7] "3L-stage_fix-dt.txt" u ($1*cyc[algsRest[i]]):(column(algsRest[i]+1)/2**(L/2)/($1*cyc[algsRest[i]])) w l dt i t word(schemes, algsRest[i])

set logscale xy
set key bottom left
set xlabel "$q/h$"
set ylabel "error"

set output "2-stage_fix-t_opt.tex"
plot [2:][:3] for [i=1:7] "2-stage_fix-t.txt" u 1:(column(algsOpt[i]+1)/2**(L/2)) w l dt i t word(schemes, algsOpt[i])

set output "3L-stage_fix-t_opt.tex"
plot [2:][:3] for [i=1:7] "3L-stage_fix-t.txt" u 1:(column(algsOpt[i]+1)/2**(L/2)) w l dt i t word(schemes, algsOpt[i])
