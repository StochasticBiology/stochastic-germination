reset

vpp(beta, delta, lambda, Lambda, nu) = lambda*(beta**2*lambda**2*Lambda**2 + delta**2*nu*(delta+nu) + beta*delta*lambda*Lambda*(delta+2*(lambda*Lambda+nu)))/(nu*(beta*lambda*Lambda + delta*nu)*(beta*lambda*Lambda + delta*(delta+nu)))
mpp(lambda, nu) = lambda/nu

set xrange [0.00001:10]
set yrange [0.1:10]
set view map
set isosamples 50
set xlabel "β"
set logscale x
set logscale y
set ylabel "Λ"
splot sqrt(vpp(x, 1., 10., y, 0.1))/mpp(10., 0.1) w pm3d
