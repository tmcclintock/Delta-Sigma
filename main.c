/*
  This is the wrapper file that
  will call the DeltaSigma code.
 */
#include <stdio.h>
#include <stdlib.h>

int read_file(FILE *fp,int N, double *data);

int main(int argc, char **argv){
  FILE *k_fp, *P_fp;
  char * line = NULL;
  size_t len = 0;
  int i,N = 0;
  int read;

  k_fp = fopen("test_data/matter_power_lin/k_h.txt","r");
  while ((read = getline(&line,&len,k_fp)) != -1){
    N++;
  }
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*k = (double*)malloc((N-1)*sizeof(double));
  read_file(k_fp,N-1,k);

  P_fp = fopen("test_data/matter_power_lin/p_k.txt","r");
  getline(&line,&len,P_fp); //header line read off
  double*P = (double*)malloc((N-1)*sizeof(double));
  read_file(P_fp,N-1,P);
  
  free(k),free(P);
  fclose(k_fp),fclose(P_fp);
}

int read_file(FILE *fp,int N, double *data){
  int i;
  for(i = 0; i < N; i++)
    fscanf(fp,"%lf",&data[i]);
  return 0;
}

