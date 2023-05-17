set style fill solid

set size ratio -1
unset xtics
unset ytics
unset border
set terminal png size 1024,1024
set output 'conf_final.png'
plot 'conf_final.dat' title "" with ellipses