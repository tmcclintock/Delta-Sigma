#include "wrapper.h"

int interface(double*k_lin,double*P_lin,int Nk_lin,
	      double*k,double*P,int Nk,
	      int NR,double Rmin,double Rmax,
	      cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs){
  int i;

  double*R=outputs->R;
  double*xi_1halo=outputs->xi_1halo;
  double*xi_mm=outputs->xi_mm;
  double*xi_lin=outputs->xi_lin;
  double*xi_2halo=outputs->xi_2halo;
  double*xi_hm=outputs->xi_hm;
  double*sigma_r=outputs->sigma_r;
  double*delta_sigma=outputs->delta_sigma;
  double*Rbins=outputs->Rbins;
  double*ave_delta_sigma=outputs->ave_delta_sigma;
  double*bias=outputs->bias;
  double*nu=outputs->nu;
  double*sigma_mis=outputs->sigma_mis;
  double*delta_sigma_mis=outputs->delta_sigma_mis;
  double*miscentered_sigma_r=outputs->miscentered_sigma_r;
  double*miscentered_delta_sigma=outputs->miscentered_delta_sigma;
  double*ave_miscentered_delta_sigma=outputs->ave_miscentered_delta_sigma;
  double*ave_delta_sigma_mis=outputs->ave_delta_sigma_mis;

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
  /*
    These two variables are for the hankel transformation
    They can be changed by hand, depending on how wacky the
    the power spectrum looks. It is difficult to test
    what combination works best for what scenarios.
  */
  int Nevals = 200;
  double h = 0.005;
  calc_xi_mm(R,NR,k,P,Nk,xi_mm,err,Nevals,h);
  if (timing){
    printf("xi_mm time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_xi_mm(R,NR,k_lin,P_lin,Nk_lin,xi_lin,err,Nevals,h);
  if (timing){
    printf("xi_lin time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_tinker_bias(&Mass,1,k_lin,P_lin,Nk_lin,bias,nu,delta,cosmo);
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
    calc_ave_delta_sigma(R,NR,delta_sigma,Nbins,R_bin_min,R_bin_max,Rbins,ave_delta_sigma);
    if (timing){
      printf("ave_delta_sigma time = %f\n",omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
  }


  if(miscentering){
    calc_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma_r,NR,sigma_mis,err,cosmo);
    if (timing){
      printf("sigma_mis time = %f\n",
	     omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
    calc_delta_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma_r,sigma_mis,
			 NR,delta_sigma_mis,err,cosmo);
    if (timing){
      printf("delta_sigma_mis time = %f\n",
	     omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
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
    if (averaging){
      calc_ave_miscentered_delta_sigma(R,NR,miscentered_delta_sigma,Nbins,R_bin_min,R_bin_max,Rbins,ave_miscentered_delta_sigma);
      if (timing){
	printf("ave_miscentered_delta_sigma time = %f\n",
	       omp_get_wtime()-time);fflush(stdout);
	time=omp_get_wtime();
      }
      /*NOTE: the ave_delta_sigma_mis is calculated with
	an identical function call as ave_miscentered_delta_sigma.
       */
      calc_ave_miscentered_delta_sigma(R,NR,delta_sigma_mis,Nbins,R_bin_min,R_bin_max,Rbins,ave_delta_sigma_mis);
      if (timing){
	printf("ave_delta_sigma_mis time = %f\n",
	       omp_get_wtime()-time);fflush(stdout);
	time=omp_get_wtime();
      }

    }
  }

  free(err);
  return 0;
}

int python_interface(double*k_lin,double*P_lin,int Nk_lin,
		     double*k_nl,double*P_nl,int Nk_nl,
		     int NR,double Rmin,double Rmax,
		     double h,double om,double ode,double ok,
		     double Mass, double concentration,
		     double Rmis, double fmis, int delta,
		     int*flow_control,int timing, int miscentering,
		     int averaging, int Nbins,
		     double R_bin_min, double R_bin_max,
		     double*R,double*xi_1halo,double*xi_mm,double*xi_lin,
		     double*xi_2halo,double*xi_hm,double*sigma_r,
		     double*delta_sigma,double*Rbins,
		     double*ave_delta_sigma,double*bias,
		     double*nu,double*sigma_mis,
		     double*delta_sigma_mis,
		     double*miscentered_sigma_r,
		     double*miscentered_delta_sigma,
		     double*ave_miscentered_delta_sigma,
		     double*ave_delta_sigma_mis){

  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->h = h;
  cosmo->om = om;
  cosmo->ode = ode;
  cosmo->ok = ok;

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->fmis=fmis;
  params->timing=timing; //1 is true
  params->miscentering=miscentering;//1 is true
  params->averaging=averaging; //1 is true
  params->Nbins=Nbins;
  params->R_bin_min=R_bin_min;
  params->R_bin_max=R_bin_max;

  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  outputs->R=R;
  outputs->xi_1halo=xi_1halo;
  outputs->xi_mm=xi_mm;
  outputs->xi_lin=xi_lin;
  outputs->xi_2halo=xi_2halo;
  outputs->xi_hm=xi_hm;
  outputs->sigma_r=sigma_r;
  outputs->delta_sigma=delta_sigma;
  outputs->bias=bias;
  outputs->nu=nu;
  outputs->sigma_mis=sigma_mis;
  outputs->delta_sigma_mis=delta_sigma_mis;
  outputs->miscentered_sigma_r=miscentered_sigma_r;
  outputs->miscentered_delta_sigma=miscentered_delta_sigma;
  outputs->Rbins=Rbins;
  outputs->ave_delta_sigma=ave_delta_sigma;
  outputs->ave_miscentered_delta_sigma=ave_miscentered_delta_sigma;
  outputs->ave_delta_sigma_mis=ave_delta_sigma_mis;

  interface(k_lin,P_lin,Nk_lin,k_nl,P_nl,Nk_nl,
	    NR,Rmin,Rmax,*cosmo,params,outputs);

  return 0;
}
		     
