#include "miscentered_ave_delta_sigma.h"

int calc_miscentered_ave_delta_sigma(double*R,int NR,
				     double*miscentered_delta_sigma,
				     int Nbins,
				     double R_bin_min,double R_bin_max,
				     double*Rbins,
				     double*miscentered_ave_delta_sigma){
  int i, status=0;
  double lRmin=log(R_bin_min),lRmax=log(R_bin_max);
  double dlR=(lRmax-lRmin)/Nbins;
  double lRlow,lRhigh,Rlow,Rhigh; //high and low bounds of current bin
#pragma omp parallel shared(R,Rbins,miscentered_delta_sigma,miscentered_ave_delta_sigma) private(i,lRlow,lRhigh,Rlow,Rhigh)
#pragma omp for
  for( i = 0; i < Nbins; i++){
    lRlow = lRmin+i*dlR;
    lRhigh = lRmin+(i+1)*dlR;
    Rlow = exp(lRlow), Rhigh = exp(lRhigh);
    Rbins[i] = (2./3.)*(Rhigh*Rhigh*Rhigh-Rlow*Rlow*Rlow)/(Rhigh*Rhigh-Rlow*Rlow); //radial mean of the bin
    //Rbins[i] = (Rhigh+Rlow)/2.; //raw mean of the bin
    status |= calc_miscentered_ave_delta_sigma_in_bin(R,NR,miscentered_delta_sigma,lRlow,lRhigh,&miscentered_ave_delta_sigma[i]);
  }
  return status;
}
