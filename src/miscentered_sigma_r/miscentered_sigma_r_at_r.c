#include "miscentered_sigma_r_at_r.h"

#define TOL1 1e-3
#define TOL2 1e-4
#define EPS 0.00001 //error offset for angular integral
#define workspace_size 8000
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

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

static double P_mis(double Rc,double Rmis_sq){
  return Rc/(Rmis_sq)*exp(-Rc*Rc/(2.0*Rmis_sq));
}//2D gaussian of width Rmis and mean 0

static int do_integral(double*sigmar_r,double*err,integrand_params*params);

static double integrand_outer(double theta,void*params);
static double integrand_inner(double lRc,void*params);

static double sigma_r_1halo_analytic(double R,double Mass,double concentration,double om,double H0,int delta);

int calc_miscentered_sigma_r_at_r(double Rp,double Mass,double concentration,
				  int delta,double Rmis,double*R,double*sigma_r,
				  int NR,double*mis_sigma_r,double*err,
				  cosmology cosmo){

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,sigma_r,NR);
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

  do_integral(mis_sigma_r,err,params);
  *mis_sigma_r *= 1./PI;
  *err *= 1./PI;
  //Factor of PI from the angular integral

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

int do_integral(double*mis_sigma_r,double*err,integrand_params*params){
  int status = 0;

  gsl_integration_workspace*workspace=params->workspace;

  gsl_function F;
  F.function=&integrand_outer;
  F.params=params;

  double result,abserr;

  status = gsl_integration_qag(&F,0,PI,TOL1,TOL1/10.,workspace_size,6,workspace,&result,&abserr);
  
  *mis_sigma_r = result; 
  *err = abserr;

  return status;
}

double integrand_outer(double theta,void*params){

  integrand_params*pars = (integrand_params*)params;
  double cos_theta = cos(theta);
  pars->cos_theta = cos_theta;

  double lrmin = pars->lrmin,lrmax = pars->lrmax;
  
  gsl_integration_workspace*workspace=pars->workspace2;
  gsl_function F;
  F.function = &integrand_inner;
  F.params = pars;

  double result,abserr;
  int status = 0;
  status = gsl_integration_qag(&F,lrmin-10,lrmax,TOL1,TOL1/10.,workspace_size,6,workspace,&result,&abserr);

  return result;
}

double integrand_inner(double lRc,void*params){
  double Rc = exp(lRc);

  integrand_params*pars = (integrand_params*)params;

  double rmin = pars->rmin,rmax = pars->rmax;
  double Rp = pars->rperp;
  double cos_theta = pars->cos_theta;
  double arg = sqrt(Rp*Rp+Rc*Rc-2*Rp*Rc*cos_theta);
  if (arg < rmin){
    double Mass = pars->Mass;
    double concentration = pars->concentration;
    int delta = pars->delta;
    cosmology cosmo = pars->cosmo;
    double om = cosmo.om;
    double h = cosmo.h;
    return Rc*P_mis(Rc,pars->Rmis_sq)*sigma_r_1halo_analytic(arg,Mass,concentration,om,h*100.,delta);
  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    return Rc*P_mis(Rc,pars->Rmis_sq)*gsl_spline_eval(spline,arg,acc);
  }else{
    return 0;
  }
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
