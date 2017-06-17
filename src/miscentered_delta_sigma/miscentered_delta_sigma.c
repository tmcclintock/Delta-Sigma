#include "miscentered_delta_sigma.h"

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
  calc_miscentered_sigma_at_r(Rp,Mass,concentration,delta,Rmis,R,sigma,NR,&answer,&error,cosmo);
  return Rp*Rp*answer;
}

int calc_miscentered_delta_sigma(double*Rp,double Mass,double concentration,
				 int delta,double Rmis,double*R,
				 double*sigma,double*miscentered_sigma,
				 int NR,double*miscentered_delta_sigma,
				 double*err,cosmology cosmo){
  int i, status = 0;
  double lrmin = log(Rp[0]);
  double inner_result=0,abserr1=0;

  if(R[0]/Rmis > 0.05){
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
    
    gsl_integration_workspace*workspace
      = gsl_integration_workspace_alloc(workspace_size);
    gsl_function F;
    F.params=params;
    F.function=&integrand_inner;
    status |= gsl_integration_qag(&F,lrmin-10,lrmin,TOL,TOL/10.,workspace_size,6,workspace,&inner_result,&abserr1);
  }else{
    //If Rmin ~< Rmis then we can use a power law approximationd
    double alpha = (log(miscentered_sigma[0])-log(miscentered_sigma[1]))/(log(Rp[0])-log(Rp[1]));
    double A = miscentered_sigma[0]/pow(Rp[0],alpha);
    inner_result = A/(alpha+2.0)*(pow(Rp[0],alpha+2.0)-pow(exp(lrmin-10),alpha+2.0));
  }
  
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_delta_sigma_at_r(Rp[i],Mass,concentration,delta,
						Rmis,R,sigma,
						miscentered_sigma,inner_result,
						NR,&miscentered_delta_sigma[i],
						&err[i],cosmo);
  }
  return status;
}
