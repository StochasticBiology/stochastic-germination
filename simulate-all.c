// Stochastic Germination calculations

// simple stochastic simulation of the ABA system under different parameterisations
// we don't even use the Gillespie algorithm here; for simplicity we just use constant time interval simulation and Poisson events based on random number checking. separate investigation shows that this doesn't introduce noteworthy issues

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()
#define NOUT 5
#define NEXPT 10
#define MAXT 100

void Simulate(int model, double beta, double Lambda)
{
  int aba, e, d;
  double dt = 0.0001;
  double betas = 0.1, betad = 0.2, delta = 1, lambda = 10, nu = 0.1, Lambdas = 2, Lambdad = 2.;
  double meanaba, meane, meand, sdaba, sde, sdd;
  int *abarec, *erec, *drec;
  FILE *fp, *fp1;
  char str[100];
  int expt;
  int i;
  double t;
  double maba, vaba, vpp, ratio;
  
  // allocate memory to record stats of individual time series
  abarec = (int*)malloc(sizeof(int)*NEXPT*MAXT);
  erec = (int*)malloc(sizeof(int)*NEXPT*MAXT);
  drec = (int*)malloc(sizeof(int)*NEXPT*MAXT);

  // the set of reactions we are considering are:
  // aba -> aba + e
  // aba -> aba + d
  // e -> e + aba
  // d -> d - aba
  // e -> 0
  // d -> 0
  // aba -> 0
  // 0 -> aba

  printf("Simulating model %i with beta = %f Lambda = %f\n", model, beta, Lambda);
  betad = betas = beta; Lambdad = Lambdas = Lambda;

  // other parameters of the model
  delta = 1; lambda = 10; nu = 0.1;
  sprintf(str, "simulate-%i-%.2f-%.2f", model, beta, Lambda);
  fp = fopen(str, "w");

  // loop through many instances of stochastic simulation
  for(expt = 0; expt < NEXPT; expt++)
    {
      // output a certain number of individual traces to files
      if(expt < NOUT)
	{
	  sprintf(str, "simulate-trace-%i-%i-%.2f-%.2f.txt", expt, model, beta, Lambda);
	  fp1 = fopen(str, "w");
	}

      // initial conditions
      aba = 10; e = 0; d = 0;

      // loop through time
      for(t = 0; t < MAXT; t += dt)
	{
	  // different models here correspond to different regulatory architectures as in the paper
	  // for each, use a simple random number check to see if an instance of this Poisson process occurs in this timestep
	  // this suffers from the flaws that the Gillespie algorithm fixes (e.g. what if multiple instances should occur in this timestep; but here we've set dt low enough and checked with other work (and theory) that these problems are negligible
	  if(model == 0)
	    {
	      if(RND < dt*betas*aba) e++;
	      if(RND < dt*betad*aba) d++;
	      if(RND < dt*delta*e) e--;
	      if(RND < dt*delta*d) d--;
	      if(RND < dt*nu*aba*(Lambdad*d+1)) aba--;
	      if(RND < dt*lambda*(Lambdas*e+1)) aba++;
	    }
	  if(model == 1)
	    {
	      if(RND < dt*nu*aba*(Lambdad*aba+1)) aba--;
	      if(RND < dt*lambda*(Lambdas*aba+1)) aba++;
	    }
	  if(model == 2)
	    {
	      if(RND < dt*betas*aba) e++;
	      if(RND < dt*delta*e) e--;
	      if(RND < dt*nu*aba*(Lambdad*e+1)) aba--;
	      if(RND < dt*lambda*(Lambdas*e+1)) aba++;
	    }

	  // this shouldn't happen, but exists as a sanity check
	  if(aba < 0) aba = 0;

	  // record state of system if we've crossed an integer timestep
	  if((int)(t) != (int)(t-dt))
	    {
	      abarec[expt*MAXT+(int)(t)] = aba;
	      erec[expt*MAXT+(int)(t)] = e;
	      drec[expt*MAXT+(int)(t)] = d;
	      if(expt < NOUT)
		fprintf(fp1, "%f %i %i %i\n", t, aba, e, d);
	    }
	}
      if(expt < NOUT)
	fclose(fp1);
    }

  // simulations are complete, so we calculate statistics across the set of simulated traces
  for(t = 0; t < MAXT; t++)
    {
      meanaba = meane = meand = 0;
      // mean levels of each species
      for(i = 0; i < NEXPT; i++)
	{
	  meanaba += abarec[i*MAXT+(int)(t)];
	  meane += erec[i*MAXT+(int)(t)];
	  meand += drec[i*MAXT+(int)(t)];
	}
      meanaba /= NEXPT; meane /= NEXPT; meand /= NEXPT;
      sdaba = sde = sdd = 0;
      // standard deviations
      for(i = 0; i < NEXPT; i++)
	{
	  sdaba += (meanaba-abarec[i*MAXT+(int)(t)])*(meanaba-abarec[i*MAXT+(int)(t)]);
	  sde += (meane-erec[i*MAXT+(int)(t)])*(meane-erec[i*MAXT+(int)(t)]);
	  sdd += (meand-drec[i*MAXT+(int)(t)])*(meand-drec[i*MAXT+(int)(t)]);
	}
      sdaba = sqrt(sdaba/(NEXPT-1)); sde = sqrt(sde/(NEXPT-1)); sdd = sqrt(sdd/(NEXPT-1));

   
      // predicted ABA statistics from theory
      maba = (betas*lambda*Lambdas - nu + sqrt(4.*betad*lambda*Lambdad*nu + (nu - betas*lambda*Lambdas)*(nu - betas*lambda*Lambdas)))/(2.*betad*Lambdad*nu);
      vaba = (lambda*(-1. + betas*lambda*Lambdas - nu) + maba*(betas*lambda*Lambdas*(-1.+(-2.+betas)*lambda*Lambdas) - (1. + 2.*betad*lambda*Lambdad)*nu - nu*nu - betad*Lambdad*nu*maba*(1.+betas*lambda*Lambdas + 3.*nu + 2.*(1.+betad)*Lambdad*nu*maba)))/(2.*(betas*lambda*Lambdas - nu - 2.*betad*Lambdad*nu*maba)*(1.+nu+betad*Lambdad*nu*maba));
      if(Lambdad == 0) { maba = lambda/nu; vaba = lambda/nu;}
      vpp = lambda*(beta*beta*lambda*lambda*Lambda*Lambda + delta*delta*nu*(delta+nu) + beta*delta*lambda*Lambda*(delta + 2*(lambda*Lambda + nu))) / (nu*(beta*lambda*Lambda + delta*nu)*(beta*lambda*Lambda + delta*(delta+nu)));
      ratio = lambda / nu;
		  
      // output to file
      fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", t, meanaba, meane, meand, sdaba, sde, sdd, vpp, ratio);
      //fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f\n", t, betas, betad, Lambdas, Lambdad, meanaba, meane, meand, sdaba, sde, sdd, maba, vaba);
    }
  fclose(fp);

  free(abarec);
  free(erec);
  free(drec);
}

int main(void)
{
  int model;
  double beta, Lambda;

  // first for Fig 2, use basic model (0) and sweep through parameters      
  model = 0;
  for(beta = 0.0; beta <= 10; beta *= 10)
    {
      for(Lambda = 0.0; Lambda <= 10; Lambda *= 10)
	{
	  Simulate(0, beta, Lambda);
	  // this line sets up the parameter loop structure after the first zero iteration
	  if(Lambda == 0) Lambda = 0.001;
	}
      // this line sets up the parameter loop structure after the first zero iteration
      if(beta == 0) beta = 0.001;
    }	

  // now for Fig 6
  // loop through different model structures
  for(model = 0; model <= 2; model++)
    {
      // loop through parameterisation, here varying betas = betad and Lambdas = Lambdad
      beta = 0.1;
      for(Lambda = 0; Lambda <= 2; Lambda += 0.5)
	{
	  Simulate(model, beta, Lambda);
	}
    }

  return 0;
}
