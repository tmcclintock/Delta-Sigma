#include "wrapper.h"

int interface(double*k,double*P,int Nk,int NR,double Rmin,double Rmax,
	      cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs){
  int i;

  double*R=outputs->R;
  double*xi_1halo=outputs->xi_1halo;
  double*xi_mm=outputs->xi_mm;
  double*xi_2halo=outputs->xi_2halo;
  double*xi_hm=outputs->xi_hm;
  double*sigma_r=outputs->sigma_r;
  double*delta_sigma=outputs->delta_sigma;
  double*bias=outputs->bias;
  double*nu=outputs->nu;
  double*miscentered_sigma_r=outputs->miscentered_sigma_r;
  double*miscentered_delta_sigma=outputs->miscentered_delta_sigma;
  
  //Used to hold integration errors for some routines
  double*err=(double*)malloc(NR*sizeof(double));

  double Mass=params->Mass;
  double concentration=params->concentration;
  double Rmis=params->Rmis;
  double fmis=params->fmis;
  int delta=params->delta;

  int Nbins=params->Nbins;
  double R_bin_min=params->R_bin_min;
  double R_bin_max=params->R_bin_max;

  int timing=params->timing;
  int miscentering=params->miscentering;
  int averaging=params->averaging;
  int*flow_control=params->flow_control;

  for(i = 0; i < NR; i++){
    double dlR = (log(Rmax) - log(Rmin))/(float)(NR-1.);
    R[i] = exp(log(Rmin) + i*dlR);
  }

  double time=omp_get_wtime();
  calc_xi_nfw(R,NR,Mass,concentration,delta,xi_1halo,cosmo);
  if (timing){
    printf("xi_nfw time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_xi_mm(R,NR,k,P,Nk,xi_mm,err);
  if (timing){
    printf("xi_mm time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_tinker_bias(&Mass,1,k,P,Nk,bias,nu,delta,cosmo);
  if (timing){
    printf("tinker_bias time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_xi_2halo(NR,xi_mm,*bias,xi_2halo);
  if (timing){
    printf("xi_2halo time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_xi_hm(NR,Mass,xi_1halo,xi_2halo,xi_hm);
  if (timing){
    printf("xi_hm time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_sigma_r(R,Mass,concentration,delta,R,xi_hm,NR,sigma_r,err,cosmo);
  if (timing){
    printf("sigma_r time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_delta_sigma(R,Mass,concentration,delta,R,sigma_r,NR,delta_sigma,err,cosmo);
  if (timing){
    printf("delta_sigma time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }

  if(averaging){
    //calc_ave_delta_sigma(R,NR,delta_sigma,Nbins,R_bin_min,R_bin_max,ave_delta_sigma,cosmo);
    if (timing){
      printf("ave_delta_sigma time = %f\n",omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
  }


  if(miscentering){
    calc_miscentered_sigma_r(R,Mass,concentration,delta,Rmis,R,sigma_r,NR,
			     miscentered_sigma_r,err,cosmo);
    if (timing){
      printf("miscentered_sigma_r time = %f\n",
	     omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
    calc_miscentered_delta_sigma(R,Mass,concentration,delta,Rmis,R,sigma_r,
				 miscentered_sigma_r,NR,miscentered_delta_sigma,
				 err,cosmo);
    if (timing){
      printf("miscentered_delta_sigma time = %f\n",
	     omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
  }


  free(err);
  return 0;
}
