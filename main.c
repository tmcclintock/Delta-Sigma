/*
  This is the wrapper file that
  will call the DeltaSigma code.
 */
#include <stdio.h>
#include <stdlib.h>
#include "src/xi_mm/xi_mm.h"

int read_file(FILE *fp,int N, double *data);

int main(int argc, char **argv){
  FILE *k_fp, *P_fp;
  char * line = NULL;
  size_t len = 0;
  int i,N = 0;
  int read;

  k_fp = fopen("test_data/matter_power_nl/k_h.txt","r");
  while ((read = getline(&line,&len,k_fp)) != -1){
    N++;
  }
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*k = (double*)malloc((N-1)*sizeof(double));
  read_file(k_fp,N-1,k);

  P_fp = fopen("test_data/matter_power_nl/p_k.txt","r");
  getline(&line,&len,P_fp); //header line read off
  double*P = (double*)malloc((N-1)*sizeof(double));
  read_file(P_fp,N-1,P);

  int NR = 100;
  double R[NR];
  double xi[NR];
  double err[NR];
  for(i = 0; i < NR; i++){
    double dlR = (log(150) - log(0.01))/(float)(NR-1.);
    R[i] = exp(log(0.01) + i*dlR);
  }
  calc_xi_mm(R,NR,k,P,N-1,xi,err);
  for(i = 0; i < NR; i++){
    printf("%e %e %e\n",R[i],xi[i],err[i]);
  }

  FILE *xiout = fopen("output/xi_mm.txt","w");
  FILE *Rout = fopen("output/R.txt","w");
  for(i = 0; i < NR; i++){
    fprintf(xiout,"%e %e\n",xi[i],err[i]);
    fprintf(Rout,"%e\n",R[i]);
  }

  free(k),free(P);
  fclose(k_fp),fclose(P_fp);
}

int read_file(FILE *fp,int N, double *data){
  int i;
  for(i = 0; i < N; i++)
    fscanf(fp,"%lf",&data[i]);
  return 0;
}

