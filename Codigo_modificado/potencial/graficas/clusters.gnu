set style fill solid

set size ratio -1

set terminal png size 1024,1024
set output 'clusters.png'
plot 'despues.dat' title "" with ellipses