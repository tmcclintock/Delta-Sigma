#include "miscentered_sigma_r_at_r.h"

#define TOL1 1e-3
#define TOL2 1e-4
#define EPS 0.001 //error offset for angular integral
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
  double Rc; //Integration variable
  double Rp2plusRc2;//Temp variable
  double RpRctimes2;//Temp variable for math optimization
}integrand_params;

static double P_mc(double Rc,double Rmis_sq){
  return Rc/(Rmis_sq)*exp(-Rc*Rc/(2.0*Rmis_sq));
}//2D gaussian of width Rmis and mean 0

static int do_integral(double*sigmar_r,double*err,integrand_params*params);

static double integrand_outer(double Rc,void*params);
static double integrand_inner_boundaries(double cos_theta,void*params);
static double integrand_inner_small_scales(double cos_theta,void*params);
static double integrand_inner_spline_scales(double cos_theta,void*params);

static double sigma_r_1halo_analytic(double R,double Mass,double concentration,double om,double H0,int delta);

int calc_miscentered_sigma_r_at_r(double Rp,double Mass,double concentration,
				  int delta,double Rmis,double*R,double*sigma_r,
				  int NR,double*mis_sigma_r,double*err,
				  cosmology cosmo){

  double h = cosmo.h, om = cosmo.om;
  double H0 = h*100.;
  double rhom = om*3.*(H0*Mpcperkm*H0*Mpcperkm)/(8.*PI*G)
    /(h*h*1e12);//SM h^2/pc^2/Mpc
  //The integral is over R, which is Mpc/h

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
  *mis_sigma_r *= rhom*2;
  *err *= rhom*2;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

int do_integral(double*mis_sigma_r,double*err,integrand_params*params){
  int status = 0;

  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  double Rp = params->rperp;
  double rmax = exp(lrmax);
  gsl_integration_workspace*workspace=params->workspace;

  gsl_function F;
  F.function=&integrand_outer;
  F.params=params;

  double result,abserr;
  status = gsl_integration_qag(&F,lrmin-10,lrmax,//log(sqrt(rmax*rmax-Rp*Rp)),
			       TOL1,TOL1/10.,workspace_size,6,workspace,&result,&abserr);

  *mis_sigma_r = result*2; //Factor of 2*PI multiplied on from 2D gaussian integral
  //, the pi is dropped from algebra
  *err = abserr;

  return status;
}

double integrand_outer(double lRc,void*params){

  integrand_params*pars = (integrand_params*)params;
  double Rc = exp(lRc);
  pars->Rc=Rc;
  double Rp = pars->rperp;
  double rmin = pars->rmin,rmax = pars->rmax;
  
  gsl_integration_workspace*workspace=pars->workspace2;
  gsl_function F;
  F.params = pars;

  double result,abserr;
  int status = 0;
  double Rp2plusRc2 = Rp*Rp+Rc*Rc;
  pars->Rp2plusRc2 = Rp2plusRc2;
  double RpRctimes2 = 2*Rp*Rc;
  pars->RpRctimes2 = RpRctimes2;
  if ((Rp2plusRc2 + RpRctimes2) < rmin*rmin){
    //Everything is below rmin
    F.function = &integrand_inner_small_scales;
  }else if((Rp2plusRc2 - RpRctimes2) > rmin*rmin && (Rp2plusRc2 + RpRctimes2) < rmax*rmax){
    //Everything is on the spline
    F.function = &integrand_inner_spline_scales;
  }else if((Rp2plusRc2 - RpRctimes2) > rmax*rmax){ 
    //F.function = &integrand_inner_large_scales;
    return 0; //Everything is past rmax (end of the spline)
  }else{
    F.function = &integrand_inner_boundaries;
  }
  //F.function = &integrand_inner_boundaries;
  status |= gsl_integration_qag(&F,0,PI,TOL2,TOL2/10.,
				//status |= gsl_integration_qag(&F,-1+EPS,1-EPS,TOL2,TOL2/10.,
				workspace_size,6,workspace,&result,&abserr);
  return Rc*P_mc(Rc,pars->Rmis_sq)*result;
}

//This is the function that needs to be broken into three parts
//Basically this will amount to putting the if statements in the outer integrand.
double integrand_inner_boundaries(double cos_theta,void*params){

  integrand_params*pars = (integrand_params*)params;

  double Rp = pars->rperp;
  double Rc = pars->Rc;
  double rmin = pars->rmin,rmax = pars->rmax;
  double Rp2plusRc2 = pars->Rp2plusRc2;
  double RpRctimes2 = pars->RpRctimes2;
  //double arg = sqrt(Rp2plusRc2 + RpRctimes2*cos_theta);
  //double arg = sqrt(pars->Rp2plusRc2 + pars->RpRctimes2*cos(cos_theta));
  double arg = sqrt(Rp2plusRc2 + RpRctimes2*cos(cos_theta));

  if (arg<rmin){
    double Mass = pars->Mass;
    double concentration = pars->concentration;
    int delta = pars->delta;
    cosmology cosmo = pars->cosmo;
    double om = cosmo.om;
    double h = cosmo.h;
    //return sigma_r_1halo_analytic(arg,Mass,concentration,om,h*100.,delta)/sqrt(1-cos_theta*cos_theta);
    return sigma_r_1halo_analytic(arg,Mass,concentration,om,h*100.,delta);

  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    //return gsl_spline_eval(spline,arg,acc)/sqrt(1-cos_theta*cos_theta);
    return gsl_spline_eval(spline,arg,acc);
  }else{ //arg > rmax
     return 0;
  }
}

double integrand_inner_small_scales(double cos_theta,void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rp = pars->rperp;
  double Rc = pars->Rc;
  double Rp2plusRc2 = pars->Rp2plusRc2;
  double RpRctimes2 = pars->RpRctimes2;
  //double arg = sqrt(Rp*Rp + Rc*Rc + 2*Rp*Rc*cos_theta);
  double arg = sqrt(Rp2plusRc2 + RpRctimes2*cos(cos_theta));

  double Mass = pars->Mass;
  double concentration = pars->concentration;
  int delta = pars->delta;
  double om = pars->cosmo.om;
  double h = pars->cosmo.h;
  return sigma_r_1halo_analytic(arg,Mass,concentration,om,h*100.,delta);
}

double integrand_inner_spline_scales(double cos_theta,void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rp = pars->rperp;
  double Rc = pars->Rc;
  double Rp2plusRc2 = pars->Rp2plusRc2;
  double RpRctimes2 = pars->RpRctimes2;
  //double arg = sqrt(Rp*Rp + Rc*Rc + 2*Rp*Rc*cos_theta);
  double arg = sqrt(Rp2plusRc2 + RpRctimes2*cos(cos_theta));

  gsl_spline*spline = pars->spline;
  gsl_interp_accel*acc = pars->acc;
  return gsl_spline_eval(spline,arg,acc);
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
