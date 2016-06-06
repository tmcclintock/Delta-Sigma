#include "xi_nfw.h"

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds
#define delta_c 1.686 //Critical collapse density
#define PI 3.141592653589793 //apple flavored

int calc_xi_nfw(double*R,int NR,double Mass,double concentration,int delta,double*xi_nfw,cosmology cosmo){
  int i;

  double h = cosmo.h, om = cosmo.om;
  double H0 = h*100.;
  double rhom = om*3.*(H0*Mpcperkm*H0*Mpcperkm)/(8.*PI*G)
    /(h*h);//SM h^2/Mpc^3
  double Rdelta = pow(Mass/(4./3.*PI*rhom*delta),1./3.);
  double Rscale = Rdelta/concentration;
  double fc = log(1.+concentration)-concentration/(1.+concentration);
  
  for(i = 0; i < NR; i++){
    xi_nfw[i] = Mass/(4.*PI*Rscale*Rscale*Rscale*fc)/(R[i]/Rscale*(1+R[i]/Rscale)*(1+R[i]/Rscale))/rhom - 1.0;
  }
  return 0;
}
