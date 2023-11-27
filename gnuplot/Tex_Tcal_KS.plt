reset

set terminal png
set output 'images/Tex_Tcal_KS.png'

set palette defined (-12 "blue", 12 "red")

set cbrange [-12:12]


set pm3d map

set title "Heatmap de Tex-Tcal"

set xlabel "x"
set ylabel "y"


splot 'dat/Tex_Tcal_KS.dat' using 1:2:3 with pm3d
