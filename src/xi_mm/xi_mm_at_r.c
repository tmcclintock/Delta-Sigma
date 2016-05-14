#include "xi_mm_at_r.h"

#define TOL 1e-8
#define workspace_size 8000
#define PI 3.141592653589793

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating xi at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace * workspace;
  double r;
  double lkmin;
  double lkmax;
}integrand_params;

static int do_integral(double*xi,double*err,integrand_params*params);

static double integrand(double lk, void*params);

int calc_xi_mm_at_r(double R,
		    double*k,double*P,
		    int N,
		    double*xi,
		    double*err){

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,N);
  gsl_spline_init(spline,k,P,N);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  
  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->r=R;
  params->lkmin=log(k[0]);
  params->lkmax=log(k[N-1]);
  
  do_integral(xi,err,params);

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int do_integral(double*xi,double*err,integrand_params*params){
  gsl_function F;
  F.function=&integrand;
  F.params=params;

  double lkmin = params->lkmin;
  double lkmax = params->lkmax;
  gsl_integration_workspace*workspace=params->workspace;

  /*
    Goal: rewrite this integral so that it goes period by period.
    This way it terminates when enough periods have completed,
    and it doesn't have to integrate over the whole P(k).
   */

  double result,abserr;
  int status;
  if(1==0){
  
    int i=1;
    double R = params->r;
    double lk0 = log(PI/R); //The k where sin(x) is 0 first.
    while(lk0 < lkmin){
      i++;
      lk0 = log(i*PI/R);
    }// Now i is the first root that is greater than lkmin

    status = gsl_integration_qag(&F,lkmin,lk0,TOL,TOL/10.,workspace_size,6,workspace,&result,&abserr);//Result now has the integral from lkmin to lk0
    double next_result=0, lk1=log((i+1)*PI/R);
    while(lk1 < lkmax && fabs((result-next_result)/result) > TOL){
      status |= gsl_integration_qag(&F,lk0,lk1,TOL,TOL/10.,workspace_size,6,workspace,&next_result,&abserr);
      //if (i<10){printf("%d %e %e %e %e %e %e\n",i,lkmin,lkmax,lk0,lk1,result,next_result);}
      fflush(stdout);
      result+=next_result;
      i++;
      lk0 = log(i*PI/R);
      lk1 = log((i+1)*PI/R);
      //if (i<10){printf("%d %e %e %e %e %e %e\n",i,lkmin,lkmax,lk0,lk1,result,next_result);}
    }
  }else{
    status = gsl_integration_qag(&F,lkmin,lkmax,TOL,TOL/10.,workspace_size,6,workspace,&result,&abserr);
  }
  
  *xi = result/(2.*PI*PI);
  *err= abserr/(2.*PI*PI);
  return status;
}

double integrand(double lk, void*params){
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double R = pars.r;
  
  double k = exp(lk);
  double x  = k*R;
  double P = gsl_spline_eval(spline,k,acc);
  return k*k*k*P*sin(x)/x;
}
