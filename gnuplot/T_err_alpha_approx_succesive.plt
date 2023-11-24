reset 

set terminal png

unset key
set output 'images/T_err_alpha_approx_succesive.png'

# Style du tracé, ici courbe bleu
set style line 1 lc rgb '#0074D9' lt 1 lw 2

# Titre du graphique
set title "f3-Uex(.,H) en norme L²"



# Étiquettes des axes
set xlabel "alpha"

set ylabel "err"

# On limite selon y pour ne pas avoir un graphe trop étendu
set yrange [0:100]  

# Tracé de la courbe 
plot 'dat/err_alpha_approx_succesive.dat' using 1:2 with lines ls 1
