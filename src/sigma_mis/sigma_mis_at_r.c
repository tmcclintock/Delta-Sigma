#include "sigma_mis_at_r.h"

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating xi at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  gsl_integration_workspace*workspace2;
  double rperp;
  double rperp_sq;
  double lrmin;
  double lrmax;
  double rmin;
  double rmax;
  cosmology cosmo;
  double Mass;
  double concentration;
  int delta;
  double Rmis; //Miscentering length
  double Rmis_sq; //Rmis^2
  double cos_theta; //Temp variable for cos(theta)
}integrand_params;

static int do_integral(double*sigmar,double*err,integrand_params*params);

static double integrand(double theta,void*params);

static double sigma_1halo_analytic(double R,double Mass,double concentration,double om,int delta);

int calc_sigma_mis_at_r(double Rp,double Mass,double concentration,
				  int delta,double Rmis,double*R,double*sigma,
				  int NR,double*mis_sigma,double*err,
				  cosmology cosmo){

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace * workspace2
    = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->workspace2=workspace2;
  params->rperp=Rp;
  params->rperp_sq=Rp*Rp;
  params->lrmin=log(R[0]);
  params->lrmax=log(R[NR-1]);
  params->rmin=R[0];
  params->rmax=R[NR-1];
  params->cosmo=cosmo;
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->Rmis_sq=Rmis*Rmis;

  do_integral(mis_sigma,err,params);
  *mis_sigma *= invPI;
  *err *= invPI;
  //Factor of PI from the angular integral. 
  //The 2 is cancelled because our integral goes from 0 to PI.

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

int do_integral(double*mis_sigma,double*err,integrand_params*params){
  int status = 0;

  //Create the workspace and the function with parameters
  gsl_integration_workspace*workspace=params->workspace;
  gsl_function F;
  F.function=&integrand;
  F.params=params;

  double result,abserr;

  status = gsl_integration_qag(&F,0,PI,MISCENTERED_TOL,MISCENTERED_TOL2,workspace_size,6,workspace,&result,&abserr);
  
  *mis_sigma = result; 
  *err = abserr;

  return status;
}

double integrand(double theta,void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rp   = pars->rperp;
  double Rmis = pars->Rmis;
  double arg  = sqrt(Rp*Rp+Rmis*Rmis-2*Rp*Rmis*cos(theta));

  double rmin = pars->rmin,rmax = pars->rmax;
  if (arg < rmin){
    double Mass = pars->Mass;
    double concentration = pars->concentration;
    int delta = pars->delta;
    cosmology cosmo = pars->cosmo;
    double om = cosmo.om;
    return sigma_1halo_analytic(arg,Mass,concentration,om,delta);
  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    return gsl_spline_eval(spline,arg,acc);
  }

  //If arg > rmax then return 0
  return 0;
}

/*
  The analytic form of Sigma(R) for an NFW profile.
 */
double sigma_1halo_analytic(double R,double Mass,double concentration,
			      double om,int delta){
  double c = concentration;
  double rhom = om*rhomconst; //SM h^2/Mpc^3
  double deltac = (delta/3.)*c*c*c/(log(1.+c)-c/(1.+c));
  double rdelta = pow(Mass/(4./3.*PI*rhom*delta),1./3.);//Mpc/h
  double rscale = rdelta/c;
  double x = R/rscale;
  double gx = 0;
  if(x<1){
    gx = (1 - 2./sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x))))/(x*x-1);
  }else if(x>1){
    gx = (1 - 2./sqrt(x*x-1)* atan(sqrt((x-1)/(1+x))))/(x*x-1);
  }
  else{ // x is approximately equal to 1
    gx = 1./3.;
  }
  return 2*rscale*deltac*rhom*gx/(1.e12);//SM h/pc^2
}
