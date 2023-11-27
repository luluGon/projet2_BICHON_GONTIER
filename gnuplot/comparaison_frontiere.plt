reset 

set terminal png

unset key
set output 'images/comparaison_frontiere.png'

# Style du tracé
set style line 1 lc  rgb "blue" lt 1 lw 1
set style line 2 lc  rgb "red" lt 2 lw 1

# Titre du graphique
set title "Comparaison frontiere "

# Étiquettes des axes
set xlabel "x"

set ylabel "y"

# On limite selon y pour ne pas avoir un graphe trop étendu
set yrange [-2:2]  

# Tracé de la courbe 
plot "dat/frontiere_alpha_opt.dat" using 1:2 with lines ls 1, "dat/Gamma_ex.dat" using 1:2 with points ls 2
