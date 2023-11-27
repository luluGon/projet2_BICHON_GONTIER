reset

set terminal png

set output 'images/Tcal_KS.png'

set palette defined (-12 "blue", 12 "red")

set cbrange [-12:12]


set pm3d map

set title "Heatmap de Tcal"

set xlabel "x"
set ylabel "y"


splot 'dat/Tcal(x,y)_KS.dat' using 1:2:3 with pm3d
