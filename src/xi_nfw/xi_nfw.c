#include "xi_nfw.h"

int calc_xi_nfw(double*R,int NR,double Mass,double concentration,int delta,double*xi_nfw,cosmology cosmo){
  int i;
  double om = cosmo.om;
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double Rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double Rscale = Rdelta/concentration;
  double fc = log(1.+concentration)-concentration/(1.+concentration);
  
  for(i = 0; i < NR; i++){
    xi_nfw[i] = Mass/(4.*M_PI*Rscale*Rscale*Rscale*fc)/(R[i]/Rscale*(1+R[i]/Rscale)*(1+R[i]/Rscale))/rhom - 1.0;
  }
  return 0;
}
