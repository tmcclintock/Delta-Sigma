#include "miscentered_sigma_at_r.h"

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
  double Rp2_cos_theta;
}integrand_params;

static int do_integral(double*sigmar_r,double*err,integrand_params*params);

static double integrand_outer(double theta,void*params);
static double integrand_inner(double lRc,void*params);

static double sigma_1halo_analytic(double R,double Mass,double concentration,double om,int delta);

int calc_miscentered_sigma_at_r(double Rp,double Mass,double concentration,
				int delta,double Rmis,double*R,double*sigma,
				int NR,double*mis_sigma,
				double*err,cosmology cosmo){

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
  *mis_sigma *= 1./PI;
  *err *= 1./PI;
  //Factor of PI from the angular integral

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

int do_integral(double*mis_sigma,double*err,integrand_params*params){
  int status = 0;

  gsl_integration_workspace*workspace=params->workspace;

  gsl_function F;
  F.function=&integrand_outer;
  F.params=params;

  double result,abserr;

  status = gsl_integration_qag(&F,0,PI,MISCENTERED_TOL,MISCENTERED_TOL2,workspace_size,6,workspace,&result,&abserr);
  
  *mis_sigma = result; 
  *err = abserr;

  return status;
}

double integrand_outer(double theta,void*params){

  integrand_params*pars = (integrand_params*)params;
  double cos_theta = cos(theta);
  pars->cos_theta = cos_theta;
  pars->Rp2_cos_theta = pars->rperp*cos_theta*2;

  double lrmin = pars->lrmin,lrmax = pars->lrmax;
  
  gsl_integration_workspace*workspace=pars->workspace2;
  gsl_function F;
  F.function = &integrand_inner;
  F.params = pars;

  double result,abserr;
  int status = 0;
  status = gsl_integration_qag(&F,lrmin-10,lrmax,MISCENTERED_TOL,MISCENTERED_TOL2,workspace_size,6,workspace,&result,&abserr);

  return result;
}

double integrand_inner(double lRc, void*params){
  double Rc = exp(lRc);
  double Rc_sq = Rc*Rc;

  integrand_params*pars = (integrand_params*)params;
  double rmin = pars->rmin,rmax = pars->rmax;
  double Rp = pars->rperp;
  double Rp_sq = pars->rperp_sq;
  double Rmis_sq = pars->Rmis_sq;
  double cos_theta = pars->cos_theta;
  double Rp2_cos_theta = pars->Rp2_cos_theta;
  double arg = sqrt(Rp_sq + Rc_sq - Rc*Rp2_cos_theta);
  
  double Rcsq_Rmissq = Rc_sq/Rmis_sq;

  /*Note: P_miscentering is a 2D gaussian.
    This is the Rc/Rmis^2 * exp() term.
   */
  if (arg < rmin){
    double Mass = pars->Mass;
    double conc = pars->concentration;
    cosmology cosmo = pars->cosmo;
    return Rcsq_Rmissq*exp(-0.5*Rcsq_Rmissq)*sigma_1halo_analytic(arg, Mass, conc, cosmo.om, pars->delta);
  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    return Rcsq_Rmissq*exp(-0.5*Rcsq_Rmissq)*gsl_spline_eval(spline, arg, acc);
  }else{
    return 0;
  }
}

double sigma_1halo_analytic(double R, double Mass, double concentration, double om, int delta){
  double c = concentration;
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double deltac = delta*0.333333333333*c*c*c/(log(1.+c)-c/(1.+c));
  double rdelta = pow(Mass/(1.333333333333*PI*rhom*delta), 0.333333333333);//Mpc/h
  double rscale = rdelta/c;
  double x = R/rscale;
  double gx = 0;
  if(x<1){
    gx = (1 - 2./sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x))))/(x*x-1);
  }else if(x>1){
    gx = (1 - 2./sqrt(x*x-1)* atan(sqrt((x-1)/(1+x))))/(x*x-1);
  }
  else{ // x is equal to 1
    gx = 0.333333333333;
  }
  return 2*rscale*deltac*rhom*gx*1.e-12;//SM h/pc^2
}
