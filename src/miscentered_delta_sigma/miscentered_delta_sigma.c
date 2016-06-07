#include "miscentered_delta_sigma.h"

#define TOL 1e-5
#define workspace_size 8000

typedef struct integrand_params{
  double lrmin;
  cosmology cosmo;
  double Mass;
  double concentration;
  int delta;
  double Rmis; //Miscentering length
  double*R;
  double*sigma_r;
  int NR;
}integrand_params;

static double integrand_inner(double lRp, void*params){
  double Rp = exp(lRp);
  integrand_params pars=*(integrand_params *)params;
  double Mass = pars.Mass;
  double concentration = pars.concentration;
  int delta = pars.delta;
  double Rmis = pars.Rmis;
  double*R = pars.R;
  double*sigma_r = pars.sigma_r;
  int NR = pars.NR;
  cosmology cosmo = pars.cosmo;
  double h = cosmo.h, om = cosmo.om;
  double answer = 0,error=0; //we put the answer here
  calc_miscentered_sigma_r_at_r(Rp,Mass,concentration,delta,Rmis,R,sigma_r,NR,&answer,&error,cosmo);
  return Rp*Rp*answer;
}

int calc_miscentered_delta_sigma(double*Rp,double Mass,double concentration,
				 int delta,double Rmis,double*R,
				 double*sigma_r,double*miscentered_sigma_r,
				 int NR,double*miscentered_delta_sigma,
				 double*err,cosmology cosmo){
  int i, status = 0;

  double lrmin = log(Rp[0]);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->lrmin=log(R[0]);
  params->cosmo=cosmo;
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->R=R;
  params->sigma_r=sigma_r;
  params->NR=NR;

  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.params=params;
  double inner_result=0,abserr1=0;
  F.function=&integrand_inner;
  double time=omp_get_wtime();
  status |= gsl_integration_qag(&F,lrmin-10,lrmin,TOL,TOL/10.,workspace_size,6,workspace,&inner_result,&abserr1);
  printf("inner time = %f\n",omp_get_wtime()-time);fflush(stdout);
  //inner_result contains the numerator of Sigma(<R), which is the costly
  //integral over the non-spline region

#pragma omp parallel shared(R,sigma_r,NR,miscentered_sigma_r,miscentered_delta_sigma,err,status)
#pragma omp for
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_delta_sigma_at_r(Rp[i],Mass,concentration,delta,
						Rmis,R,sigma_r,
						miscentered_sigma_r,inner_result,
						NR,&miscentered_delta_sigma[i],
						&err[i],cosmo);
  }
}
