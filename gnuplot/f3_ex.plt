reset 

set terminal png

unset key
set output 'images/f3_ex.png'

# Style du tracé, ici courbe bleu
set style line 1 lc rgb '#0074D9' lt 1 lw 2

# Titre du graphique
set title "f3 exacte"



# Étiquettes des axes
set xlabel "x"

set ylabel "f3(x)"

 

# Tracé de la courbe 
plot 'dat/f3ex.dat' using 1:2 with lines ls 1
