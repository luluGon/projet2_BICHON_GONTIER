reset 

set terminal png

unset key
set output 'images/frontiere_alpha_opt.png'

# Style du tracé, ici courbe noir
set style line 1 lc rgb 'black' lt 1 lw 2

# Titre du graphique
set title "Frontiere Gamma approché"

# Étiquettes des axes
set xlabel "x"

set ylabel "y"

# On limite selon y pour ne pas avoir un graphe trop étendu
set yrange [-2:2]  

# Tracé de la courbe 
plot "dat/frontiere_alpha_opt.dat" using 1:2 with lines ls 1
