#include "miscentered_delta_sigma_at_r.h"

#define TOL 1e-3
#define workspace_size 8000
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating mis_sigma_r at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  double lrperp;
  double lrmin;
  double lrmax;
  double rmin;
  double rmax;
  cosmology cosmo;
  double Mass;
  double concentration;
  int delta;
  double Rmis; //Miscentering length
  double*R;
  double*sigma_r;
  int NR;
}integrand_params;

static int do_integral(double*miscentered_delta_sigma,double*err,integrand_params*params);

static double integrand1(double R, void*params);
static double integrand2(double R, void*params);

static double sigma_r_1halo_analytic(double R,double Mass,double concentration,
				     double om,double H0,int delta);


int calc_miscentered_delta_sigma_at_r(double Rp,double Mass,
				      double concentration,int delta,
				      double Rmis,double*R,double*sigma_r,
				      double*miscentered_sigma_r,int NR,
				      double*miscentered_delta_sigma,
				      double*err,cosmology cosmo){

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,miscentered_sigma_r,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->lrperp=log(Rp);
  params->lrmin=log(R[0]);
  params->lrmax=log(R[NR-1]);
  params->rmin=R[0];
  params->rmax=R[NR-1];
  params->cosmo=cosmo;
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->R=R;
  params->sigma_r=sigma_r;
  params->NR=NR;

  do_integral(miscentered_delta_sigma,err,params);
  *miscentered_delta_sigma *= 2./Rp/Rp;
  *miscentered_delta_sigma -= gsl_spline_eval(spline,Rp,acc);
  *err *= 2./Rp/Rp;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int do_integral(double*miscentered_delta_sigma,double*err,integrand_params*params){
  int status = 0;

  gsl_function F;
  F.params=params;

  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  double lRp = params->lrperp;
  gsl_integration_workspace*workspace=params->workspace;

  double result1=0,abserr1=0;
  double result2=0,abserr2=0;

  F.function=&integrand1;
  status |= gsl_integration_qag(&F,lrmin-10,lrmin,TOL,TOL/10.,workspace_size,6,workspace,&result1,&abserr1);

  F.function=&integrand2;
  status |= gsl_integration_qag(&F,lrmin,lRp,TOL,TOL/10.,workspace_size,6,workspace,&result2,&abserr2);

  *miscentered_delta_sigma = result1 + result2;
  *err= abserr1 + abserr2;
  return status;
}

double integrand1(double lR, void*params){
  double Rp = exp(lR);
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
  calc_miscentered_sigma_r_at_r(Rp,Mass,concentration,delta,Rmis,R,sigma_r,NR,&answer,&error,cosmo);return Rp*Rp*answer;
  //return Rp*Rp*sigma_r_1halo_analytic(Rp,Mass,concentration,om,h*100.,delta);
}

double integrand2(double lR, void*params){
  double R = exp(lR);
  integrand_params pars=*(integrand_params *)params;
  gsl_spline*spline = pars.spline;//Sigma_R(R) spline
  gsl_interp_accel*acc = pars.acc;
  return R*R*gsl_spline_eval(spline,R,acc);
}

double sigma_r_1halo_analytic(double R,double Mass,double concentration,
			      double om,double H0,int delta){
  double c = concentration;
  double rhom = om*3.*(H0*H0*Mpcperkm*Mpcperkm)/(8.*PI*G)
    /(H0/100.*H0/100.);//SM h^2/Mpc^3
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
