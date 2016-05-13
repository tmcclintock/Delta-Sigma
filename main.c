/*
  This is the wrapper file that
  will call the DeltaSigma code.
 */
#include <stdio.h>
#include <stdlib.h>
#include "src/sigma_r/sigma_r.h"
#include "src/xi_hm/xi_hm.h"
#include "src/xi_mm/xi_mm.h"
#include "src/tinker_bias/tinker_bias.h"
#include "src/cosmology/cosmology.h"

int read_file(FILE *fp,int N, double *data);

int main(int argc, char **argv){
  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->h = 0.7;
  cosmo->om = 0.3;
  cosmo->ode = 0.7;
  cosmo->ok = 0.0;

  FILE *k_fp, *P_fp;
  char * line = NULL;
  size_t len = 0;
  int i,N = -1;
  int read;

  k_fp = fopen("test_data/matter_power_nl/k_h.txt","r");
  while ((read = getline(&line,&len,k_fp)) != -1){
    N++;
  }
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*k = (double*)malloc((N)*sizeof(double));
  read_file(k_fp,N,k);

  P_fp = fopen("test_data/matter_power_nl/p_k.txt","r");
  getline(&line,&len,P_fp); //header line read off
  double*P = (double*)malloc((N)*sizeof(double));
  read_file(P_fp,N,P);

  int NR = 300;//N;
  double R[NR];
  double xi_mm[NR];
  double err[NR];
  for(i = 0; i < NR; i++){
    double dlR = (log(200) - log(0.01))/(float)(NR-1.);
    R[i] = exp(log(0.01) + i*dlR);
  }
  calc_xi_mm(R,NR,k,P,N,xi_mm,err);

  int NM = 100;
  double M[NM];
  double bias[NM];
  double nu[NM];
  for(i = 0; i < NM; i++){
    double dlM = (16.0 - 12.0)/(float)(NM-1.);
    M[i] = pow(10,(12.0 + i*dlM));
  }
  calc_tinker_bias(M,NM,k,P,N,bias,nu,200,*cosmo);

  double Mass = M[NM/2];
  double test_bias = bias[NM/2];
  double concentration = 4.0*pow(Mass/5.e14,-0.1);//Bad M-c relation
  double xi_nfw[NR];
  calc_xi_nfw(R,NR,Mass,concentration,200,xi_nfw,*cosmo);

  double xi_hm[NR];
  calc_xi_hm(NR,Mass,xi_nfw,xi_mm,test_bias,xi_hm);

  double sigma_r[NR];
  double sigma_r_err[NR];
  calc_sigma_r(R,Mass,concentration,200,R,xi_hm,NR,sigma_r,sigma_r_err,*cosmo);

  FILE *xi_mm_out = fopen("output/xi_mm.txt","w");
  FILE *xi_nfw_out = fopen("output/xi_nfw.txt","w");
  FILE *xi_hm_out = fopen("output/xi_hm.txt","w");
  FILE *sigma_r_out = fopen("output/sigma_r.txt","w");

  FILE *Rout = fopen("output/R.txt","w");
  for(i = 0; i < NR; i++){
    fprintf(sigma_r_out,"%e\n",sigma_r[i]);
    fprintf(xi_hm_out,"%e\n",xi_hm[i]);
    fprintf(xi_nfw_out,"%e\n",xi_nfw[i]);
    fprintf(xi_mm_out,"%e %e\n",xi_mm[i],err[i]);
    fprintf(Rout,"%e\n",R[i]);
  }
  FILE *bias_out = fopen("output/bias.txt","w");
  FILE *nu_out = fopen("output/nu.txt","w");
  FILE *M_out = fopen("output/M.txt","w");
  for(i = 0; i < NM; i++){
    fprintf(bias_out,"%e\n",bias[i]);
    fprintf(nu_out,"%e\n",nu[i]);
    fprintf(M_out,"%e\n",M[i]);
  }

  free(k),free(P);
  free(cosmo);
  fclose(k_fp),fclose(P_fp);
  fclose(xi_mm_out),fclose(Rout);
  fclose(bias_out),fclose(nu_out),fclose(M_out);
}

int read_file(FILE *fp,int N, double *data){
  int i;
  for(i = 0; i < N; i++)
    fscanf(fp,"%lf",&data[i]);
  return 0;
}

