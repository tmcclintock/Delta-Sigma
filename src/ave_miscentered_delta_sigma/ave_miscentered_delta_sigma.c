#include "ave_miscentered_delta_sigma.h"

int calc_ave_miscentered_delta_sigma(double*R,int NR,
				     double*miscentered_delta_sigma,
				     int Nbins,
				     double R_bin_min,double R_bin_max,
				     double*Rbins,
				     double*ave_miscentered_delta_sigma){
  int i, status=0;
  double lRmin=log(R_bin_min),lRmax=log(R_bin_max);
  double dlR=(lRmax-lRmin)/Nbins;
  double lRlow,lRhigh,Rlow,Rhigh; //high and low bounds of current bin
  for( i = 0; i < Nbins; i++){
    lRlow = lRmin+i*dlR;
    lRhigh = lRmin+(i+1)*dlR;
    Rlow = exp(lRlow), Rhigh = exp(lRhigh);
    Rbins[i] = 0.66666666*(Rhigh*Rhigh*Rhigh-Rlow*Rlow*Rlow)/(Rhigh*Rhigh-Rlow*Rlow); //radial mean of the bin
    status |= calc_ave_miscentered_delta_sigma_in_bin(R,NR,miscentered_delta_sigma,lRlow,lRhigh,&ave_miscentered_delta_sigma[i]);
  }
  return status;
}
