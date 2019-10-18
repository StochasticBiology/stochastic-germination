# svg 800,400
reset
set multiplot
set size 0.33,1

#plot "outspan.txt" u 6:12:1 w labels tc rgbcolor "#FF0000", "" u 9:(sqrt($13)):1 w labels tc rgbcolor "#00FF00", x
set logscale
set origin 0,0
set border 3 lw 2
set grid
set xtics nomirror
set ytics nomirror
set xlabel "Stochastic simulation"
set ylabel "Analytic prediction"
set xrange [1:1e5]
set yrange [1:1e5]
set origin 0,0
set key bottom right
set xtics 1, 100, 1e5
set label 1 at graph 0.05, graph 0.95 "A." font "Arial Black, 22"
plot x notitle lc rgbcolor "#000000", "outspan.txt" u 6:12 ls 1 title "Mean", "" u 9:(sqrt($13)) pt 2 lc rgbcolor "#0000BB" title "SD", x notitle
set origin 0.33,0
unset key
set label 1 at graph 0.05, graph 0.95 "B." font "Arial Black, 22"
plot x notitle lc rgbcolor "#000000", "outspan-direct.txt" u 6:12 ls 1 title "Mean", "" u 9:(sqrt($13)) pt 2 lc rgbcolor "#0000BB" title "SD"
set origin 0.66,0
set label 1 at graph 0.05, graph 0.95 "C." font "Arial Black, 22"
plot x notitle lc rgbcolor "#000000", "outspan-single1.txt" u 6:12 ls 1 title "Mean", "" u 9:(sqrt($13)) pt 2 lc rgbcolor "#0000BB" title "SD"

