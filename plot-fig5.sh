#svg 800,400
reset
set multiplot
unset key
set size 0.5,1
set origin 0.5,0
set style fill solid
set border 3 lw 2
set xtics nomirror
set ytics nomirror
set grid
set boxwidth 0.25
set yrange [-0.02:*]
set xrange [-0.05:*]
set xtics 0.2
set xlabel "Mean germination proportion"
set ylabel "Noise in transformed\ngermination proportion"
set label 1 at graph 0.86, graph 0.95 "B." font "Arial Black, 24"
plot "fig5-data-i.txt" u 5:3:($6-$5):($4-$3) w vectors lw 2 filled head, "" u 5:3 w p ls 1 ps 1.5 pt 7, "" u ($5+$1*($6-$5)/10):($3+$1*($4-$3)/10):2 w labels 

set origin 0.,0
set boxwidth 0.25
set yrange [0:*]
set xtics 10
set ylabel "Noise in ABA level"
set xlabel "ABA (pg / seed)"
set label 1 at graph 0.86, graph 0.95 "A." font "Arial Black, 24"
plot "fig5-data-ii.txt" u 5:3:($6-$5):($4-$3) w vectors lw 2 filled head, "" u 5:3 w p ls 1 ps 1.5 pt 7, "" u ($5+($6-$5)/2):($3+($4-$3)/2):2 w labels 
