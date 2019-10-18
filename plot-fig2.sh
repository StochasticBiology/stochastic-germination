reset
set multiplot
set size 0.33,0.5
unset key

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

set origin 0,0.
set yrange [0:300]
set trange [0:200]
mean = lambda/nu
eta = etav(0,1)
set label 1 at graph 0.05, graph 0.95 "β = 0, Λ = 1"
plot "simulate-trace-0-0-0.00-1.00.txt" u 1:2 w l ls 1, "simulate-trace-1-0-0.00-1.00.txt" u 1:2 w l ls 2, "simulate-trace-2-0-0.00-1.00.txt" u 1:2 w l ls 3, "simulate-trace-3-0-0.00-1.00.txt" u 1:2 w l ls 4, "simulate-trace-4-0-0.00-1.00.txt" u 1:2 w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6
set origin 0.33,0.25
set yrange [0:300]
set trange [0:200]
mean = lambda/nu
eta = etav(0.1,1)
set label 1 at graph 0.05, graph 0.95 "β = 0.1, Λ = 1"
set parametric
plot "simulate-trace-0-0-0.10-1.00.txt" u 1:2 w l ls 1, "simulate-trace-1-0-0.10-1.00.txt" u 1:2 w l ls 2, "simulate-trace-2-0-0.10-1.00.txt" u 1:2 w l ls 3, "simulate-trace-3-0-0.10-1.00.txt" u 1:2 w l ls 4, "simulate-trace-4-0-0.10-1.00.txt" u 1:2 w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6 
set origin 0.66,0
set yrange [0:300]
set trange [0:200]
mean = lambda/nu
eta = etav(1,1)
set label 1 at graph 0.05, graph 0.95 "β = 1, Λ = 1"
plot "simulate-trace-0-0-1.00-1.00.txt" u 1:2 w l ls 1, "simulate-trace-1-0-1.00-1.00.txt" u 1:2 w l ls 2, "simulate-trace-2-0-1.00-1.00.txt" u 1:2 w l ls 3, "simulate-trace-3-0-1.00-1.00.txt" u 1:2 w l ls 4, "simulate-trace-4-0-1.00-1.00.txt" u 1:2 w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6 

##########

set origin 0,0.5
set yrange [0:300]
set trange [0:200]
mean = lambda/nu
eta = etav(0.1,0)
set label 1 at graph 0.05, graph 0.95 "β = 0.1, Λ = 0"
plot "simulate-trace-0-0-0.10-0.00.txt" u 1:2 w l ls 1, "simulate-trace-1-0-0.10-0.00.txt" u 1:2 w l ls 2, "simulate-trace-2-0-0.10-0.00.txt" u 1:2 w l ls 3, "simulate-trace-3-0-0.10-0.00.txt" u 1:2 w l ls 4, "simulate-trace-4-0-0.10-0.00.txt" u 1:2 w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6
set origin 0.66,0.5
set yrange [0:300]
set trange [0:200]
mean = lambda/nu
eta = etav(0.1,10)
set label 1 at graph 0.05, graph 0.95 "β = 0.1, Λ = 10"
plot "simulate-trace-0-0-0.10-10.00.txt" u 1:2 w l ls 1, "simulate-trace-1-0-0.10-10.00.txt" u 1:2 w l ls 2, "simulate-trace-2-0-0.10-10.00.txt" u 1:2 w l ls 3, "simulate-trace-3-0-0.10-10.00.txt" u 1:2 w l ls 4, "simulate-trace-4-0-0.10-10.00.txt" u 1:2 w l ls 5, 100-1000*gauss(t, mean, eta*mean),t ls 6 
