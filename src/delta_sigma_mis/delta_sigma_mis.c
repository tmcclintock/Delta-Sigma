#include "delta_sigma_mis.h"

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
  double*sigma;
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
  double*sigma = pars.sigma;
  int NR = pars.NR;
  cosmology cosmo = pars.cosmo;
  double h = cosmo.h, om = cosmo.om;
  double answer = 0,error=0; //we put the answer here
  calc_sigma_mis_at_r(Rp,Mass,concentration,delta,Rmis,R,sigma,NR,&answer,&error,cosmo);
  return Rp*Rp*answer;
}

int calc_delta_sigma_mis(double*Rp,double Mass,double concentration,
			 int delta,double Rmis,double*R,
			 double*sigma,double*sigma_mis,
			 int NR,double*delta_sigma_mis,
			 double*err,cosmology cosmo){
  int i, status = 0;
  double lrmin = log(Rp[0]);
  double inner_result=0,abserr1=0;

  if(R[0]/Rmis > 0.1){
    integrand_params*params=malloc(sizeof(integrand_params));
    params->lrmin=log(R[0]);
    params->cosmo=cosmo;
    params->Mass=Mass;
    params->concentration=concentration;
    params->delta=delta;
    params->Rmis=Rmis;
    params->R=R;
    params->sigma=sigma;
    params->NR=NR;
    
    gsl_integration_workspace * workspace
      = gsl_integration_workspace_alloc(workspace_size);
    gsl_function F;
    F.params=params;
    F.function=&integrand_inner;
    double time=omp_get_wtime();
    status |= gsl_integration_qag(&F,lrmin-10,lrmin,TOL,TOL/10.,workspace_size,6,workspace,&inner_result,&abserr1);
    //inner_result contains the numerator of Sigma(<R), which is the costly
    //integral over the non-spline region
  }else{
    //If Rmin ~< Rmis then we can use a power law approximation pretty well
    double alpha = (log(sigma_mis[0])-log(sigma_mis[1]))
      /(log(Rp[0])-log(Rp[1]));
    double A = sigma_mis[0]/pow(Rp[0],alpha);
    inner_result = A/(alpha+2.0)*(pow(Rp[0],alpha+2.0)-pow(exp(lrmin-10),alpha+2.0));
  }

#pragma omp parallel shared(R,sigma,NR,sigma_mis,delta_sigma_mis,err,status)
#pragma omp for
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_delta_sigma_at_r(Rp[i],Mass,concentration,delta,
						Rmis,R,sigma,
						sigma_mis,inner_result,
						NR,&delta_sigma_mis[i],
						&err[i],cosmo);
  }

  //Success
  return status;
}