#svg 800,600
reset
set multiplot
set size 0.33

lambda = 10
nu = 0.1
delta = 1
vaa(beta, Lambda) = lambda*(beta**2*lambda**2*Lambda**2 + delta**2*nu*(delta+nu) + beta*delta*lambda*Lambda*(delta + 2.*(lambda*Lambda + nu)))/(nu*(beta*lambda*Lambda + delta*nu)*(beta*lambda*Lambda+delta*(delta+nu)))
etav(beta, Lambda) = sqrt(vaa(beta,Lambda))/(lambda/nu)

gauss(x, mean, sd) = 1./sqrt(2.*3.141*sd*sd)*exp(-(x-mean)**2/(2.*sd*sd))

set style line 1 lw 2 lc rgbcolor "#FF0000"
set style line 2 lw 2 lc rgbcolor "#FF3333"
set style line 3 lw 2 lc rgbcolor "#FF6666"
set style line 4 lw 2 lc rgbcolor "#FF9999"
set style line 5 lw 2 lc rgbcolor "#FFAAAA"
set style line 6 lw 4 lc rgbcolor "#8888FF"

unset key

set xtics 20
set ytics 100
set parametric
set border 3 lw 2
set grid
set xtics nomirror
set ytics nomirror
#set logscale y
set xlabel "Time"
set ylabel "ABA level"

mean = lambda/nu
eta = etav(0.1,0)
etad = etav(0.1,0)
set yrange [0:300]
set trange [0:200]

set origin 0,0.66
set label 1 at graph 0.05,graph 0.9 "A." font "Arial Black, 24"
set label 2 at graph 0.7,graph 0.95 "Λ = 0"
plot "simulate-trace-0-0-0.10-0.00.txt" w l ls 1, "simulate-trace-1-0-0.10-0.00.txt" w l ls 2, "simulate-trace-2-0-0.10-0.00.txt" w l ls 3, "simulate-trace-3-0-0.10-0.00.txt" w l ls 4, "simulate-trace-4-0-0.10-0.00.txt" w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6
set origin 0.,0.33
unset label 2
set label 1 at graph 0.05,graph 0.9 "B." font "Arial Black, 24"
plot "simulate-trace-0-1-0.10-0.00.txt" w l ls 1, "simulate-trace-1-1-0.10-0.00.txt" w l ls 2, "simulate-trace-2-1-0.10-0.00.txt" w l ls 3, "simulate-trace-3-1-0.10-0.00.txt" w l ls 4, "simulate-trace-4-1-0.10-0.00.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6
set origin 0.,0.0
set label 1 at graph 0.05,graph 0.9 "C." font "Arial Black, 24"
plot "simulate-trace-0-2-0.10-0.00.txt" w l ls 1, "simulate-trace-1-2-0.10-0.00.txt" w l ls 2, "simulate-trace-2-2-0.10-0.00.txt" w l ls 3, "simulate-trace-3-2-0.10-0.00.txt" w l ls 4, "simulate-trace-4-2-0.10-0.00.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6

unset label 1
eta = etav(0.1,0.5)
set origin 0.33,0.66
set label 2 at graph 0.7,graph 0.95 "Λ = 0.5"
plot "simulate-trace-0-0-0.10-0.50.txt" w l ls 1, "simulate-trace-1-0-0.10-0.50.txt" w l ls 2, "simulate-trace-2-0-0.10-0.50.txt" w l ls 3, "simulate-trace-3-0-0.10-0.50.txt" w l ls 4, "simulate-trace-4-0-0.10-0.50.txt" w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6
set origin 0.33,0.33
unset label 2
plot "simulate-trace-0-1-0.10-0.50.txt" w l ls 1, "simulate-trace-1-1-0.10-0.50.txt" w l ls 2, "simulate-trace-2-1-0.10-0.50.txt" w l ls 3, "simulate-trace-3-1-0.10-0.50.txt" w l ls 4, "simulate-trace-4-1-0.10-0.50.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6
set origin 0.33,0.0
plot "simulate-trace-0-2-0.10-0.50.txt" w l ls 1, "simulate-trace-1-2-0.10-0.50.txt" w l ls 2, "simulate-trace-2-2-0.10-0.50.txt" w l ls 3, "simulate-trace-3-2-0.10-0.50.txt" w l ls 4, "simulate-trace-4-2-0.10-0.50.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6

eta = etav(0.1,2.0)
set origin 0.66,0.66
set label 2 at graph 0.7,graph 0.95 "Λ = 2"
plot "simulate-trace-0-0-0.10-2.00.txt" w l ls 1, "simulate-trace-1-0-0.10-2.00.txt" w l ls 2, "simulate-trace-2-0-0.10-2.00.txt" w l ls 3, "simulate-trace-3-0-0.10-2.00.txt" w l ls 4, "simulate-trace-4-0-0.10-2.00.txt" w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6
set origin 0.66,0.33
unset label 2
plot "simulate-trace-0-1-0.10-2.00.txt" w l ls 1, "simulate-trace-1-1-0.10-2.00.txt" w l ls 2, "simulate-trace-2-1-0.10-2.00.txt" w l ls 3, "simulate-trace-3-1-0.10-2.00.txt" w l ls 4, "simulate-trace-4-1-0.10-2.00.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6
set origin 0.66,0.0
plot "simulate-trace-0-2-0.10-2.00.txt" w l ls 1, "simulate-trace-1-2-0.10-2.00.txt" w l ls 2, "simulate-trace-2-2-0.10-2.00.txt" w l ls 3, "simulate-trace-3-2-0.10-2.00.txt" w l ls 4, "simulate-trace-4-2-0.10-2.00.txt" w l ls 5, 100-1000*gauss(t, mean, etad*mean),t ls 6
