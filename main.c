/*
This is a main function that tests the wrapper::python_interface() function,
in order to figure out if there are memory leaks.
 */

#include<stdio.h>
#include<stdlib.h>
#include"src/wrapper/wrapper.h"

int main(){
  printf("test main\n");
  //TODO: allocate everything and pass it in to the function

  //Step 1: read in the k and Ps

  //Step 2: define the model parameters

  //Step 3: define the radial values and the binning
  int NR = 200;
  int Nbins = 15;
  double*R=(double*)malloc(NR*sizeof(double));
  double*xi1h=(double*)malloc(NR*sizeof(double));
  double*ximm=(double*)malloc(NR*sizeof(double));
  double*xilin=(double*)malloc(NR*sizeof(double));
  double*xi2h=(double*)malloc(NR*sizeof(double));
  double*xihm=(double*)malloc(NR*sizeof(double));
  double*sigma=(double*)malloc(NR*sizeof(double));
  double*deltasigma=(double*)malloc(NR*sizeof(double));
  double*ads=(double*)malloc(Nbins*sizeof(double));
  double*bias=(double*)malloc(sizeof(double));
  double*nu=(double*)malloc(sizeof(double));
  double*sigma_mis=(double*)malloc(NR*sizeof(double));
  double*deltasigma_mis=(double*)malloc(NR*sizeof(double));
  double*mis_sigma=(double*)malloc(NR*sizeof(double));
  double*mis_deltasigma=(double*)malloc(NR*sizeof(double));
  double*ads_mis=(double*)malloc(Nbins*sizeof(double));
  double*mis_ads=(double*)malloc(Nbins*sizeof(double));

  //Step 4: allocate arrays for the IO

  //Step 5: call the python interface function

  //Step 6: free everthing allocated in step 4
  free(R);
  free(xi1h);
  free(ximm);
  free(xilin);
  free(xi2h);
  free(xihm);
  free(sigma);
  free(deltasigma);
  free(ads);
  free(bias);
  free(nu);
  free(sigma_mis);
  free(deltasigma_mis);
  free(mis_sigma);
  free(mis_deltasigma);
  free(ads_mis);
  free(mis_ads);

  return 0;
}
