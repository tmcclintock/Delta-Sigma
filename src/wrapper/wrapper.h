/*
  This is an interface between the DeltaSigma code and 
  a piece of code that can pass in a cosmology and 
  a power spectrum.
  In addition, the next thing to be implemented will be 
  a boolean array that provides flow control with regards
  to which parts of the code to run.
*/

#include "../cosmology/cosmology.h"

#ifndef INTERFACE
#define INTERFACE

typedef struct interface_parameters{
  double Mass;
  double concentration;
  double lnc;
  double fmis;
  int*flow_control;
}interface_parameters;
#endif

int interface(double*k,double*P,int N);
//,cosmology cosmo, interface_parameters params);
