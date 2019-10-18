reset
unset key

set multiplot
mpp(lambda, nu, Lambdas, Lambdad) = (lambda*Lambdas-nu+sqrt(4*lambda*Lambdad*nu+(nu-lambda*Lambdas)**2))/(2.*Lambdad*nu)
vpp(lambda, nu, Lambdas, Lambdad) = -(lambda+lambda*Lambdas*mpp(lambda, nu, Lambdas, Lambdad) + nu*mpp(lambda, nu, Lambdas, Lambdad) + Lambdad*nu*mpp(lambda, nu, Lambdas, Lambdad)**2)/(2.*(lambda*Lambdas-nu-2.*Lambdad*nu*mpp(lambda, nu, Lambdas, Lambdad)))
eta(lambda, nu, Lambdas, Lambdad) = sqrt(vpp(lambda, nu, Lambdas, Lambdad))/mpp(lambda, nu, Lambdas, Lambdad)
scaledeta(lambda, nu, Lambdas, Lambdad) = sqrt(mpp(lambda, nu, Lambdas, Lambdad))*eta(lambda, nu, Lambdas, Lambdad)
loggedmpp(lambda, nu, Lambdas, Lambdad) = log10(mpp(lambda, nu, Lambdas, Lambdad))

mpp1(lambda, nu, beta, Lambdas, Lambdad) = (beta*lambda*Lambdas - nu + sqrt((beta*lambda*Lambdas - nu)**2 + 4.*beta*lambda*Lambdad*nu))/(2.*beta*Lambdad*nu)
vpp1(lambda, nu, beta, Lambdas, Lambdad) = (beta*lambda**2*Lambdas*(1.+mpp1(lambda, nu, beta, Lambdas, Lambdad)*(beta-2)*Lambdas) - mpp1(lambda, nu, beta, Lambdas, Lambdad)*nu*(1.+nu+mpp1(lambda, nu, beta, Lambdas, Lambdad)*beta*Lambdad*(1+(3+2*mpp1(lambda, nu, beta, Lambdas, Lambdad)*(beta+1)*Lambdad)*nu)) - lambda*(1+nu+mpp1(lambda, nu, beta, Lambdas, Lambdad)*beta*(Lambdas+Lambdad*(2+mpp1(lambda, nu, beta, Lambdas, Lambdad)*(beta-4)*Lambdas)*nu)))/(2.*(beta*lambda*Lambdas-nu-2.*mpp1(lambda, nu, beta, Lambdas, Lambdad)*beta*Lambdad*nu)*(1.+nu+mpp1(lambda, nu, beta, Lambdas, Lambdad)*beta*Lambdad*nu))
eta1(lambda, nu, beta, Lambdas, Lambdad) = sqrt(vpp1(lambda, nu, beta, Lambdas, Lambdad))/mpp1(lambda, nu, beta, Lambdas, Lambdad)
scaledeta1(lambda, nu, beta, Lambdas, Lambdad) = sqrt(mpp1(lambda, nu, beta, Lambdas, Lambdad))*eta1(lambda, nu, beta, Lambdas, Lambdad)
loggedmpp1(lambda, nu, beta, Lambdas, Lambdad) = log10(mpp1(lambda, nu, beta, Lambdas, Lambdad))

tlambda = 10
tbeta = 0.1
tnu = 0.1
tLambdas = 0.1
tLambdad = 0.1

set logscale x
set logscale y
set xrange [0.001:100]
set yrange [0.001:100]
set view map
set isosamples 30

set cbrange [0:0.7]

set size 0.33,0.3

##############

set origin 0.,0.75
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot eta(tlambda, tnu, x, y) w pm3d

set origin 0.,0.5
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λs"
splot eta1(tlambda, tnu, x, y, tLambdad) w pm3d


set origin 0.,0.0
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot eta1(tlambda, tnu, tbeta, x, y) w pm3d

set origin 0.,0.25
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λd"
splot eta1(tlambda, tnu, x, tLambdas, y) w pm3d

############

set cbrange [1:5]
set cbtics 1

set origin 0.33,0.75
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot scaledeta(tlambda, tnu, x, y) w pm3d

set origin 0.33,0.5
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λs"
splot scaledeta1(tlambda, tnu, x, y, tLambdad) w pm3d

set origin 0.33,0.25
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λd"
splot scaledeta1(tlambda, tnu, x, tLambdas, y) w pm3d

set origin 0.33,0.0
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot scaledeta1(tlambda, tnu, tbeta, x, y) w pm3d

###########

set cbrange [0:5]
set palette defined (0 "#FF0000", 2 "#FFFFFF", 5 "#0000FF")

set origin 0.66,0.0
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot loggedmpp1(tlambda, tnu, tbeta, x, y) w pm3d


set origin 0.66,0.25
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λd"
splot loggedmpp1(tlambda, tnu, x, tLambdas, y) w pm3d

set origin 0.66,0.75
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot loggedmpp(tlambda, tnu, x, y) w pm3d


set origin 0.66,0.5
set xlabel offset 0,0.5 "β"
set ylabel offset 2 "Λs"
splot loggedmpp1(tlambda, tnu, x, y, tLambdad) w pm3d

##

