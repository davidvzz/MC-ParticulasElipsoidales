set style fill solid
set xrange [0:1]
set yrange [0:1]
set size ratio -1

set terminal png size 1024,1024
set output 'clusters.png'
plot 'despues.dat' title "" with ellipses