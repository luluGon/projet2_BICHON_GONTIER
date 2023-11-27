reset 

set terminal png

set key
set output 'images/f3_comparaison.png'


set style line 1 lc  rgb "blue" lt 1 lw 2
set style line 2 lc  rgb "red" lt 1 lw 2
set style line 3 lc  rgb "green" lt 1 lw 2
set style line 4 lc  rgb "brown" lt 1 lw 2

# Étiquettes des axes
set xlabel "x"

set ylabel "f3(x)"

 

# Tracé de la courbe 
plot 'dat/f3_adov1.dat' using 1:2 with lines ls 1 title "adomainv1",'dat/f3ex.dat' using 1:2 with lines ls 2 title "f3exacte",'dat/f3_app_succ.dat' using 1:2 with lines ls 3 title "app_succ",'dat/f3_ker_sep.dat' using 1:2 with lines ls 4 title "ker_sep"

