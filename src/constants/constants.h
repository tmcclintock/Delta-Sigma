//Libraries that every file should have
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//First define tolerances for integrations
#ifndef TOLERANCES
#define TOLERANCES
#define TOL 1e-5
#define AVE_TOL 1e-8 //Used for taking averages
#define BIAS_TOL 1e-8 //Used for the tinker bias
#define MISCENTERED_TOL 1e-2 // Used for miscentering
#define MISCENTERED_TOL2 MISCENTERED_TOL*0.1
#endif

//Now define physical constants
#ifndef PHYSICAL_CONSTANTS
#define PHYSICAL_CONSTANTS
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds
#define delta_c 1.686 //Critical collapse density
#define rhomconst 2.775808e+11//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#endif

//Now define constants
#ifndef CONSTANTS
#define CONSTANTS
#define PI 3.141592653589793
#define invPI 0.318309886183790 // 1/PI
#define workspace_size 8000
#endif
