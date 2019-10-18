# svg 800,1000
reset
set multiplot
set size 0.33,0.22
unset key

maba(betas, Lambdas, betad, Lambdad, lambda, nu) = (betas*lambda*Lambdas - nu + sqrt(4.*betad*lambda*Lambdad*nu + (nu - betas*lambda*Lambdas)*(nu - betas*lambda*Lambdas)))/(2.*betad*Lambdad*nu)
vaba(betas, Lambdas, betad, Lambdad, lambda, nu) = (lambda*(-1. + betas*lambda*Lambdas - nu) + maba(betas, Lambdas, betad, Lambdad, lambda, nu)*(betas*lambda*Lambdas*(-1.+(-2.+betas)*lambda*Lambdas) - (1. + 2.*betad*lambda*Lambdad)*nu - nu*nu - betad*Lambdad*nu*maba(betas, Lambdas, betad, Lambdad, lambda, nu)*(1.+betas*lambda*Lambdas + 3.*nu + 2.*(1.+betad)*Lambdad*nu*maba(betas, Lambdas, betad, Lambdad, lambda, nu))))/(2.*(betas*lambda*Lambdas - nu - 2.*betad*Lambdad*nu*maba(betas, Lambdas, betad, Lambdad, lambda, nu))*(1.+nu+betad*Lambdad*nu*maba(betas, Lambdas, betad, Lambdad, lambda, nu)))
eta(betas, Lambdas, betad, Lambdad, lambda, nu) = sqrt(vaba(betas, Lambdas, betad, Lambdad, lambda, nu))/maba(betas, Lambdas, betad, Lambdad, lambda, nu)
scaledeta(betas, Lambdas, betad, Lambdad, lambda, nu) = (sqrt(vaba(betas, Lambdas, betad, Lambdad, lambda, nu))/maba(betas, Lambdas, betad, Lambdad, lambda, nu))*sqrt(maba(betas, Lambdas, betad, Lambdad, lambda, nu))
loggedmaba(betas, Lambdas, betad, Lambdad, lambda, nu) = log10(maba(betas, Lambdas, betad, Lambdad, lambda, nu))

y1 = 0
y2 = 0.166
y3 = 0.333
y4 = 0.499
y5 = 0.666
y6 = 0.833

set view map

set cbtics 0.2
set xrange [0.001:10]
#set yrange [0.1:10]
set yrange [0.01:100]
set isosamples 30
set logscale y
tlambda = 10
tnu = 0.1
tbetas = 0.1
tLambdas = 1
tbetad = 0.1
tLambdad = 1

######### top 2



############# lower 4

set origin 0,y6
set logscale x
set cbrange [0:0.7]
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot eta(tbetas, x, tbetad, y, tlambda, tnu) w pm3d
set origin 0,y5
set logscale x
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "βd"
splot eta(x, tLambdas, y, tLambdad, tlambda, tnu) w pm3d
set origin 0,y4
set logscale x
set cbrange [0:0.7]
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λs"
splot eta(x, y, tbetad, tLambdad, tlambda, tnu) w pm3d
set origin 0,y3
set logscale x
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λd"
splot eta(tbetas, tLambdas, x, y, tlambda, tnu) w pm3d
set origin 0,y2
set logscale x
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λd"
splot eta(x, tLambdas, tbetad, y, tlambda, tnu) w pm3d
set origin 0,y1
set logscale x
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λs"
splot eta(tbetas, y, x, tLambdad, tlambda, tnu) w pm3d

#set cbrange [0:0.5]
#set cbrange [0:10]
#set palette defined (0 "#FF0000", 1 "#FFFFFF", 10 "#0000FF")
set cbtics 1
set cbrange [1:5]

#set cbrange [0:0.5]
#set cbrange [0:10]
#set palette defined (0 "#FF0000", 1 "#FFFFFF", 10 "#0000FF")
set cbrange [1:5]
set origin 0.33,y6
set logscale x
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot scaledeta(tbetas, x, tbetad, y, tlambda, tnu) w pm3d
set origin 0.33,y5
set logscale x
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "βd"
splot scaledeta(x, tLambdas, y, tLambdad, tlambda, tnu) w pm3d
set origin 0.33,y4
set logscale x
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λs"
splot scaledeta(x, y, tbetad, tLambdad, tlambda, tnu) w pm3d
set origin 0.33,y3
set logscale x
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λd"
splot scaledeta(tbetas, tLambdas, x, y, tlambda, tnu) w pm3d
set origin 0.33,y2
set logscale x
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λd"
splot scaledeta(x, tLambdas, tbetad, y, tlambda, tnu) w pm3d
set origin 0.33,y1
set logscale x
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λs"
splot scaledeta(tbetas, y, x, tLambdad, tlambda, tnu) w pm3d


set cbrange [0:5]
set palette defined (0 "#FF0000", 2 "#FFFFFF", 5 "#0000FF")
set origin 0.66,y6
set xlabel offset 0,0.5 "Λs"
set ylabel offset 2 "Λd"
splot loggedmaba(tbetas, x, tbetad, y, tlambda, tnu) w pm3d
set origin 0.66,y5
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "βd"
splot loggedmaba(x, tLambdas, y, tLambdad, tlambda, tnu) w pm3d
set origin 0.66,y4
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λs"
splot loggedmaba(x, y, tbetad, tLambdad, tlambda, tnu) w pm3d
set origin 0.66,y3
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λd"
splot loggedmaba(tbetas, tLambdas, x, y, tlambda, tnu) w pm3d
set origin 0.66,y2
set xlabel offset 0,0.5 "βs"
set ylabel offset 2 "Λd"
splot loggedmaba(x, tLambdas, tbetad, y, tlambda, tnu) w pm3d
set origin 0.66,y1
set xlabel offset 0,0.5 "βd"
set ylabel offset 2 "Λs"
splot loggedmaba(tbetas, y, x, tLambdad, tlambda, tnu) w pm3d
