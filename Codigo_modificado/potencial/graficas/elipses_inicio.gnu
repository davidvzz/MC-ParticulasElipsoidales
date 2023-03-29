set style fill solid
set xrange [0:1]
set yrange [0:1]
set size ratio -1
unset xtics
unset ytics
unset border
set terminal png size 1024,1024
set output 'conf_inicial.png'
plot 'conf_inicial.dat' title "" with ellipses