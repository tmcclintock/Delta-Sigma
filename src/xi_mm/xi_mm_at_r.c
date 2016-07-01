#include "xi_mm_at_r.h"

#define TOL 1e-8  //integral tolerance
#define TOL2 1e-4 //periodicity tolerance; Xi is right to 0.01%
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
  double kmin = exp(lkmin), kmax = exp(lkmax);
  gsl_integration_workspace*workspace=params->workspace;

  /*
    Goal: rewrite this integral so that it goes period by period.
    This way it terminates when enough periods have completed,
    and it doesn't have to integrate over the whole P(k).
   */

  double result=0,abserr,next_result=1e99;
  int status;
  if(1==1){
  
    int i=1;
    double R = params->r;
    double lPR = log(PI/R);
    double lk0 = lPR;
    while(lk0 < lkmin){
      i+=1;
      lk0 = log(i)+lPR;
    }// Now i is the first root that is greater than lkmin
    if (lk0 > lkmax){
      lk0 = lkmax;
    }

    double k_period = exp(lk0); //period of the kernal in k-space
    double periods = (kmax-kmin)/k_period;
    int step = (int)(log10(periods)+1)*2;
    //printf("%e %e %e %e %d\n",kmin,kmax,k_period,periods,step);

    status = gsl_integration_qag(&F,lkmin,lk0,TOL,TOL/10.,workspace_size,6,workspace,&result,&abserr);
    double lk1=log(i+step)+lPR;//log((i+step)*PI/R);
    while(lk1 < lkmax && fabs(next_result/result) > TOL2){
      status |= gsl_integration_qag(&F,lk0,lk1,TOL,TOL/10.,workspace_size,6,workspace,&next_result,&abserr);
      result+=next_result;
      i+=step;
      lk0 = lk1;
      lk1 = log(i+step)+lPR;
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
