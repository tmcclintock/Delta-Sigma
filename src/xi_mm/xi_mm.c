#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979

//The power spectrum
double get_P(double x,double R,double*k,double*P,int Nk,gsl_spline*Pspl,gsl_interp_accel*acc){
  double ki = x/R;
  double kmin = k[0];
  double kmax = k[Nk-1];
  double alpha,A;
  if (ki < kmin){
    alpha = log(P[1]/P[0])/log(k[1]/k[0]);
    A = P[0]/pow(k[0],alpha);
    return A*pow(ki,alpha);
  }else if (ki > kmax){
    alpha = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
    A = P[Nk-1]/pow(k[Nk-1],alpha);
    return A*pow(ki,alpha);
  }// Assume power laws at ends
  return gsl_spline_eval(Pspl,ki,acc);
}

double calc_corr_at_R(double R,double*k,double*P,int Nk,int N,double h){
  double zero,psi,x,t,dpsi,f,PIsinht;
  double PI_h = PI/h;
  double PI_2 = PI/2.;
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(Pspl,k,P,Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    zero = i+1;
    psi = h*zero*tanh(sinh(h*zero)*PI_2);
    x = psi*PI_h;
    t = h*zero;
    PIsinht = PI*sinh(t);
    dpsi = (PI*t*cosh(t)+sinh(PIsinht))/(1+cosh(PIsinht));
    if (dpsi!=dpsi) dpsi=1.0;
    f = x*get_P(x,R,k,P,Nk,Pspl,acc);
    sum += f*sin(x)*dpsi;
  }

  gsl_spline_free(Pspl),gsl_interp_accel_free(acc);
  return sum/(R*R*R*PI*2);
}

int calc_xi_mm(double*R,int NR,double*k,double*P,int Nk,double*xi,double*err,int N, double h){
  int i;

  //#pragma omp parallel for shared(R,NR,k,P,Nk,xi,err,N,h) private(i)
  for(i=0;i<NR;i++)
    xi[i] = calc_corr_at_R(R[i],k,P,Nk,N,h);
  
  return 0;
}
