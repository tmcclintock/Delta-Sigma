/*
This is a main function that tests the wrapper::python_interface() function,
in order to figure out if there are memory leaks.
 */

#include<stdio.h>
#include<stdlib.h>
#include"src/wrapper/wrapper.h"

int read_file(FILE *fp,int N, double *data);

int main(){
  printf("Starting the test main.\n");
  
  //Step 1: read in the k and Ps
  FILE*k_fp, *P_fp;
  char*line = NULL;
  size_t len = 0;
  int Nk = -1, Nklin = -1;
  int read;
  //Get the nl data
  k_fp = fopen("test_data/matter_power_nl/k_h.txt","r");
  while ((read = getline(&line, &len, k_fp)) != -1){Nk++;}
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*knl = (double*)malloc(Nk*sizeof(double));
  read_file(k_fp, Nk, knl);
  fclose(k_fp);
  P_fp = fopen("test_data/matter_power_nl/p_k.txt","r");
  getline(&line,&len,P_fp); //header line read off
  double*Pnl = (double*)malloc(Nk*sizeof(double));
  read_file(P_fp, Nk, Pnl);
  fclose(P_fp);
  //Get the lin data
  k_fp = fopen("test_data/matter_power_lin/k_h.txt","r");
  while ((read = getline(&line, &len, k_fp)) != -1){Nklin++;}
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*klin = (double*)malloc(Nklin*sizeof(double));
  read_file(k_fp, Nklin, klin);
  fclose(k_fp);
  P_fp = fopen("test_data/matter_power_lin/p_k.txt","r");
  getline(&line, &len, P_fp); //header line read off
  double*Plin = (double*)malloc(Nklin*sizeof(double));
  read_file(P_fp, Nklin, Plin);
  fclose(P_fp);
  printf("Power spectra read in.\n");

  //Step 2: define the model parameters
  double h   = 0.7; //km/s/Mpc/100
  double om  = 0.3;
  double ode = 1.0 - om;
  double ok  = 0.0;
  double Mass = 1e14; //Msun/h
  double conc = 4*pow(Mass/5e14, -0.1); //Arbitrary
  double fmis = 0.22;
  int delta = 200;
  int*flow = (int*)malloc(sizeof(int));
  int miscentering = 1;
  int averaging    = 1;
  printf("Model defined.\n");

  //Step 3: define the radial values and the binning
  int NR = 200;
  int Nbins = 15;
  double Rmin = 0.01;
  double Rmax = 200.0;
  double Rmis = 0.24;
  double R_bin_min = 0.01;
  double R_bin_max = 200.0;
  printf("Radii and bins defined.\n");

  //Step 4: allocate arrays for the IO
  double*R=(double*)malloc(NR*sizeof(double));
  double*xi1h=(double*)malloc(NR*sizeof(double));
  double*ximm=(double*)malloc(NR*sizeof(double));
  double*xilin=(double*)malloc(NR*sizeof(double));
  double*xi2h=(double*)malloc(NR*sizeof(double));
  double*xihm=(double*)malloc(NR*sizeof(double));
  double*sigma=(double*)malloc(NR*sizeof(double));
  double*deltasigma=(double*)malloc(NR*sizeof(double));
  double*Rbins=(double*)malloc(Nbins*sizeof(double));
  double*ads=(double*)malloc(Nbins*sizeof(double));
  double*bias=(double*)malloc(sizeof(double));
  double*nu=(double*)malloc(sizeof(double));
  double*sigma_mis=(double*)malloc(NR*sizeof(double));
  double*deltasigma_mis=(double*)malloc(NR*sizeof(double));
  double*mis_sigma=(double*)malloc(NR*sizeof(double));
  double*mis_deltasigma=(double*)malloc(NR*sizeof(double));
  double*ads_mis=(double*)malloc(Nbins*sizeof(double));
  double*mis_ads=(double*)malloc(Nbins*sizeof(double));
  printf("Arrays allocated.\n");

  //Step 5: call the python interface function a bunch of times
  int i;
  for(i = 0; i < 50; i++){
    python_interface(klin, Plin, Nklin,
		     knl, Pnl, Nk,
		     NR, Rmin, Rmax, 
		     h, om, ode, ok, 
		     Mass, conc,
		     Rmis, fmis, delta,
		     flow, miscentering,
		     averaging, Nbins,
		     R_bin_min, R_bin_max,
		     R, xi1h, ximm, xilin,
		     xi2h, xihm, sigma,
		     deltasigma, Rbins,
		     ads, bias,
		     nu, sigma_mis, deltasigma_mis, 
		     mis_sigma, mis_deltasigma, 
		     mis_ads, ads_mis);
  }
  printf("Python interface called.\n");

  //Step 6: free everthing allocated in step 4
  free(line);
  free(knl);
  free(klin);
  free(Pnl);
  free(Plin);
  free(R);
  free(xi1h);
  free(ximm);
  free(xilin);
  free(xi2h);
  free(xihm);
  free(sigma);
  free(deltasigma);
  free(Rbins);
  free(ads);
  free(bias);
  free(nu);
  free(sigma_mis);
  free(deltasigma_mis);
  free(mis_sigma);
  free(mis_deltasigma);
  free(ads_mis);
  free(mis_ads);
  free(flow);
  printf("Allocated values freed.\n");

  return 0;
}

int read_file(FILE *fp, int N, double *data){
  int i;
  for(i = 0; i < N; i++)
    fscanf(fp,"%lf",&data[i]);
  return 0;
}
