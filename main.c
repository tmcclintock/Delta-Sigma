/*
  This is the wrapper file that
  will call the DeltaSigma code.
 */
#include <stdio.h>
#include <stdlib.h>
#include "src/wrapper/wrapper.h"

int read_file(FILE *fp,int N, double *data);

int main(int argc, char **argv){
  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->h = 0.7;
  cosmo->om = 0.3;
  cosmo->ode = 0.7;
  cosmo->ok = 0.0;

  FILE *k_lin_fp, *P_lin_fp;
  FILE *k_nl_fp, *P_nl_fp;
  char * line = NULL;
  size_t len = 0;
  int i, N_lin = -1, N_nl = -1;
  int read;

  k_lin_fp = fopen("test_data/matter_power_lin/k_h.txt","r");
  while ((read = getline(&line,&len,k_lin_fp)) != -1){
    N_lin++;
  }
  rewind(k_lin_fp);
  read = getline(&line,&len,k_lin_fp); //header line read off
  double*k_lin = (double*)malloc((N_lin)*sizeof(double));
  read_file(k_lin_fp,N_lin,k_lin);

  P_lin_fp = fopen("test_data/matter_power_lin/p_k.txt","r");
  read = getline(&line,&len,P_lin_fp); //header line read off
  double*P_lin = (double*)malloc((N_lin)*sizeof(double));
  read_file(P_lin_fp,N_lin,P_lin);


  k_nl_fp = fopen("test_data/matter_power_nl/k_h.txt","r");
  while ((read = getline(&line,&len,k_nl_fp)) != -1){
    N_nl++;
  }
  rewind(k_nl_fp);
  read = getline(&line,&len,k_nl_fp); //header line read off
  double*k_nl = (double*)malloc((N_nl)*sizeof(double));
  read_file(k_nl_fp,N_nl,k_nl);

  P_nl_fp = fopen("test_data/matter_power_nl/p_k.txt","r");
  read = getline(&line,&len,P_nl_fp); //header line read off
  double*P_nl = (double*)malloc((N_nl)*sizeof(double));
  read_file(P_nl_fp,N_nl,P_nl);

  FILE *M_fp = fopen("test_data/mass_list.dat","r");
  int NM = 0;
  while ((read = getline(&line,&len,M_fp)) != -1){
    NM++;
  }
  rewind(M_fp);
  read = getline(&line,&len,M_fp); //header line read off
  double*M = (double*)malloc((NM)*sizeof(double));
  read_file(M_fp,NM,M);
  double*bias_arr = (double*)malloc((NM)*sizeof(double));
  double*nu_arr = (double*)malloc((NM)*sizeof(double));
  for(i = 0; i < NM; i++){
    double dlM = (16.0 - 12.0)/(float)(NM-1.);
    M[i] = pow(10,(12.0 + i*dlM));
  }
  calc_tinker_bias(M,NM,k_lin,P_lin,N_lin,bias_arr,nu_arr,200,*cosmo);

  double Mass = 1e14;
  double concentration = 4.0*pow(Mass/5.e14,-0.1);//Bad M-c relation
  double bias;
  double nu;

  int NR = 300;
  double Rmin = 0.01, Rmax = 200; //Mpc/h
  double*R=(double*)malloc(NR*sizeof(double));
  double*xi_mm=(double*)malloc(NR*sizeof(double));
  double*xi_nfw=(double*)malloc(NR*sizeof(double));
  double*xi_2halo=(double*)malloc(NR*sizeof(double));
  double*xi_hm=(double*)malloc(NR*sizeof(double));
  double*sigma_r=(double*)malloc(NR*sizeof(double));
  double*delta_sigma=(double*)malloc(NR*sizeof(double));
  double*miscentered_sigma_r=(double*)malloc(NR*sizeof(double));
  double*miscentered_delta_sigma=(double*)malloc(NR*sizeof(double));

  int Nbins = 15;
  double R_bin_min = 0.01, R_bin_max = 200; //Mpc/h
  double*Rbins=(double*)malloc(Nbins*sizeof(double));
  double*ave_delta_sigma=(double*)malloc(Nbins*sizeof(double));
  double*miscentered_ave_delta_sigma=(double*)malloc(Nbins*sizeof(double));

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=200;
  params->Rmis=0.3;
  params->fmis=0.23;
  params->timing=1; //1 is true
  params->miscentering=1;//1 is true
  params->averaging=1; //1 is true
  params->Nbins=Nbins;
  params->R_bin_min=R_bin_min;
  params->R_bin_max=R_bin_max;

  outputs->R=R;
  outputs->xi_1halo=xi_nfw;
  outputs->xi_mm=xi_mm;
  outputs->xi_2halo=xi_2halo;
  outputs->xi_hm=xi_hm;
  outputs->sigma_r=sigma_r;
  outputs->delta_sigma=delta_sigma;
  outputs->bias=&bias;
  outputs->nu=&nu;
  outputs->miscentered_sigma_r=miscentered_sigma_r;
  outputs->miscentered_delta_sigma=miscentered_delta_sigma;
  outputs->Rbins=Rbins;
  outputs->ave_delta_sigma=ave_delta_sigma;
  outputs->miscentered_ave_delta_sigma=miscentered_ave_delta_sigma;

  interface(k_lin,P_lin,N_lin,k_nl,P_nl,N_nl,
	    NR,Rmin,Rmax,*cosmo,params,outputs);

  FILE*Rout = fopen("output/R.txt","w");
  FILE*xi_mm_out = fopen("output/xi_mm.txt","w");
  FILE*xi_nfw_out = fopen("output/xi_nfw.txt","w");
  FILE*xi_2h_out = fopen("output/xi_2halo.txt","w");
  FILE*xi_hm_out = fopen("output/xi_hm.txt","w");
  FILE*sigma_r_out = fopen("output/sigma_r.txt","w");
  FILE*delta_sigma_out = fopen("output/delta_sigma.txt","w");
  FILE*miscentered_sigma_r_out;
  FILE*miscentered_delta_sigma_out;
  FILE*ave_delta_sigma_out;
  FILE*Rbins_out;
  FILE*miscentered_ave_delta_sigma_out;
  if(params->averaging){
    ave_delta_sigma_out = fopen("output/ave_delta_sigma.txt","w");
    Rbins_out = fopen("output/Rbins.txt","w");
  }
  if(params->miscentering){
    miscentered_sigma_r_out = fopen("output/miscentered_sigma_r.txt","w");
    miscentered_delta_sigma_out = fopen("output/miscentered_delta_sigma.txt","w");
    if(params->averaging)
      miscentered_ave_delta_sigma_out = fopen("output/miscentered_ave_delta_sigma.txt","w");
  }

  //Write out the averaged values
  for(i = 0; i < Nbins; i++){
    if(params->averaging){
      fprintf(ave_delta_sigma_out,"%e\n",ave_delta_sigma[i]);
      fprintf(Rbins_out,"%e\n",Rbins[i]);
      if(params->miscentering){
	fprintf(miscentered_ave_delta_sigma_out,"%e\n",
		miscentered_ave_delta_sigma[i]);
      }
    }
  }

  for(i = 0; i < NR; i++){
    if(params->miscentering){
      fprintf(miscentered_delta_sigma_out,"%e\n",miscentered_delta_sigma[i]);
      fprintf(miscentered_sigma_r_out,"%e\n",miscentered_sigma_r[i]);
    }
    fprintf(delta_sigma_out,"%e\n",delta_sigma[i]);
    fprintf(sigma_r_out,"%e\n",sigma_r[i]);
    fprintf(xi_hm_out,"%e\n",xi_hm[i]);
    fprintf(xi_2h_out,"%e\n",xi_2halo[i]);
    fprintf(xi_nfw_out,"%e\n",xi_nfw[i]);
    fprintf(xi_mm_out,"%e\n",xi_mm[i]);
    fprintf(Rout,"%e\n",R[i]);
  }
  FILE*bias_out = fopen("output/bias.txt","w");
  FILE*nu_out = fopen("output/nu.txt","w");
  FILE*M_out = fopen("output/M.txt","w");
  for(i = 0; i < NM; i++){
    fprintf(bias_out,"%e\n",bias_arr[i]);
    fprintf(nu_out,"%e\n",nu_arr[i]);
    fprintf(M_out,"%e\n",M[i]);
  }

  free(k_lin),free(P_lin);
  free(k_nl),free(P_nl);

  free(R),free(xi_mm),free(xi_nfw),free(xi_2halo);
  free(xi_hm),free(sigma_r),free(delta_sigma);
  free(cosmo);
  fclose(k_lin_fp),fclose(P_lin_fp);
  fclose(k_nl_fp),fclose(P_nl_fp);
  fclose(xi_mm_out),fclose(Rout);
  fclose(bias_out),fclose(nu_out),fclose(M_out);
  fclose(xi_2h_out),fclose(xi_nfw_out);
  fclose(xi_hm_out),fclose(sigma_r_out),fclose(delta_sigma_out);
  if(params->averaging){
    fclose(ave_delta_sigma_out);
    fclose(Rbins_out);
  }
  if(params->miscentering){
    fclose(miscentered_sigma_r_out),fclose(miscentered_delta_sigma_out);
    if(params->averaging){
      fclose(miscentered_ave_delta_sigma_out);
    }
  }
  free(params);
}

int read_file(FILE*fp,int N, double *data){
  int i, garbage;
  for(i = 0; i < N; i++)
    garbage = fscanf(fp,"%lf",&data[i]);
  return 0;
}

