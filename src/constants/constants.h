//First define tolerances for integrations
#ifndef TOLERANCES
#define TOLERANCES
#define AVE_TOL 1e-8
#define BIAS_TOL 1e-8
#define TOL 1e-5
#endif

//Now define physical constants
#ifndef PHYSICAL_CONSTANTS
#define PHYSICAL_CONSTANTS
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds
#define delta_c 1.686 //Critical collapse density
#endif

//Now define constants
#ifndef CONSTANTS
#define CONSTANTS
#define PI 3.141592653589793
#define workspace_size 8000
#endif
