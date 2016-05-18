/*
  This module will calculate the weak lensing signal DeltaSigma
  including the miscentering component.
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "src/wrapper/wrapper.h"

const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * dist = DISTANCES_SECTION;

typedef struct weak_lensing_config{
  int delta;
  double Rmin, Rmax;
  double bin_min, bin_max;
  int NR, Nbins;
  int timing=0,miscentering=0,averaging=0;
}weak_lensing_config;

void*setup(c_datablock*options){
  DATABLOCK_STATUS status=0;
  weak_lensing_config*config=malloc(sizeof(weak_lensing_config));
  status |= c_datablock_get_int(options,OPTION_SECTION,"delta",&(config->delta));
  status|=c_datablock_get_double(options,OPTION_SECTION,"Rmin",&(config->Rmin));
  status |= c_datablock_get_double(options,OPTION_SECTION,"Rmax",&(config->Rmax));
  status |= c_datablock_get_int(options,OPTION_SECTION,"NR",&(config->NR));
  status|=c_datablock_get_double(options,OPTION_SECTION,"bin_min",&(config->bin_min));
  status |= c_datablock_get_double(options,OPTION_SECTION,"bin_max",&(config->bin_max));
  status |= c_datablock_get_int(options,OPTION_SECTION,"Nbins",&(config->num_bins));
  status |= c_datablock_get_int(options,OPTION_SECTION,"timing",&(config->timing));
  status |= c_datablock_get_int(options,OPTION_SECTION,"miscentering",&(config->miscentering));
  status |= c_datablock_get_int(options,OPTION_SECTION,"averaging",&(config->averaging));
  if(status){fprintf(stderr,"Error in weak lensing setup.\n");exit(status);}
  return config;
}

int execute(c_datablock*block,void*config_in){
  int i,j,l;
  DATABLOCK_STATUS status=0;
  weak_lensing_config*config=(weak_lensing_config*)config_in;
  int delta=config->delta;
  double Rmin=config->Rmin, Rmax=config->Rmax;
  double bin_min=config->bin_min, bin_max=config->bin_max;
  int NR=config->NR, Nbins=config->Nbins;
  int timing=config->timing, miscentering=config->miscentering, averaging=config->averaging;

  //Acquire the redshift samples
  double *z;
  int NZ;
  status |= c_datablock_get_double_array_1d(block,"matter_power_lin","z",&z,&NZ);
  if(status){fprintf(stderr,"Error on reading in z.\n");exit(status);}

  //Read in the linear and nonlinear matter power spectra.
  double *k_lin;
  int NK_lin, ndims_lin;
  status |= c_datablock_get_double_array_1d(block,"matter_power_lin","k_h",&k,&NK_lin);
  if(status){fprintf(stderr,"Error on reading in k_h lin.\n");exit(status);}
  status |= c_datablock_get_array_ndim(block,"matter_power_lin","p_k",&ndims_lin);
  if(status){fprintf(stderr,"Error on reading in dimensions of PK lin.\n");exit(status);}
  int extents_lin[ndims_lin];
  double*PK_lin = (double *)malloc(NZ*NK_lin*sizeof(double));
  if (ndims_lin == 2){
    status |= c_datablock_get_double_array_shape(block,"matter_power_lin","p_k",
						 ndims_lin,extents_lin);
    if(status){fprintf(stderr,"Error on reading in 2D PK lin.\n");exit(status);}
    status |= c_datablock_get_double_array(block,"matter_power_lin",
					   "p_k",(double*)PK_lin,ndims_lin,extents_lin);
  }else{
    status |= c_datablock_get_double_array_1d(block,"matter_power_lin","p_k",&PK_lin,&NK_lin);
    if(status){fprintf(stderr,"Error on reading in 1D PK lin.\n");exit(status);}
  }
  if(status){fprintf(stderr,"Error on reading in PK lin.\n");exit(status);}
  double *k_nl;
  int NK_nl, ndims_nl;
  status |= c_datablock_get_double_array_1d(block,"matter_power_nl","k_h",&k,&NK_nl);
  if(status){fprintf(stderr,"Error on reading in k_h nl.\n");exit(status);}
  status |= c_datablock_get_array_ndim(block,"matter_power_nl","p_k",&ndims_nl);
  if(status){fprintf(stderr,"Error on reading in dimensions of PK nl.\n");exit(status);}
  int extents_nl[ndims_nl];
  double*PK_nl = (double *)malloc(NZ*NK_nl*sizeof(double));
  if (ndims_nl == 2){
    status |= c_datablock_get_double_array_shape(block,"matter_power_nl","p_k",
						 ndims_nl,extents_nl);
    if(status){fprintf(stderr,"Error on reading in 2D PK nl.\n");exit(status);}
    status |= c_datablock_get_double_array(block,"matter_power_nl",
					   "p_k",(double*)PK_nl,ndims_nl,extents_nl);
  }else{
    status |= c_datablock_get_double_array_1d(block,"matter_power_nl","p_k",&PK_nl,&NK_nl);
    if(status){fprintf(stderr,"Error on reading in 1D PK nl.\n");exit(status);}
  }
  if(status){fprintf(stderr,"Error on reading in PK nl.\n");exit(status);}

  //Create the arrays to populate and write to the block
  int extents_bias[] = {NZ};
  double*bias=(double*)malloc(NZ*sizeof(double));
  double*nu=(double*)malloc(NZ*sizeof(double));

  int ndim_out = 2;
  int extents[] = {NZ,NR};
  double**R=(double**)malloc(NZ*NR*sizeof(double*));
  double**xi_mm=(double**)malloc(NZ*NR*sizeof(double*));
  double**xi_nfw=(double**)malloc(NZ*NR*sizeof(double*));
  double**xi_2halo=(double**)malloc(NZ*NR*sizeof(double*));
  double**xi_hm=(double**)malloc(NZ*NR*sizeof(double*));
  double**sigma_r=(double**)malloc(NZ*NR*sizeof(double*));
  double**delta_sigma=(double**)malloc(NZ*NR*sizeof(double*));
  double**miscentered_sigma_r=(double**)malloc(NZ*NR*sizeof(double*));
  double**miscentered_delta_sigma=(double**)malloc(NZ*NR*sizeof(double*));
  for(i = 0; i < NZ; i++){
    R[i]=(double*)malloc(NR*sizeof(double));
    xi_mm[i]=(double*)malloc(NR*sizeof(double));
    xi_nfw[i]=(double*)malloc(NR*sizeof(double));
    xi_2halo[i]=(double*)malloc(NR*sizeof(double));
    xi_hm[i]=(double*)malloc(NR*sizeof(double));
    sigma_r[i]=(double*)malloc(NR*sizeof(double));
    delta_sigma[i]=(double*)malloc(NR*sizeof(double));
    miscentered_sigma_r[i]=(double*)malloc(NR*sizeof(double));
    miscentered_delta_sigma[i]=(double*)malloc(NR*sizeof(double));
  }

  int extents_ave[] = {NZ,Nbins};
  double**Rbins=(double**)malloc(Nbins*sizeof(double*));
  double**ave_delta_sigma=(double**)malloc(Nbins*sizeof(double*));
  double**miscentered_ave_delta_sigma=(double**)malloc(Nbins*sizeof(double*));
  for(i = 0 i < NZ; i++){
    Rbins=(double*)malloc(Nbins*sizeof(double));
    ave_delta_sigma=(double*)malloc(Nbins*sizeof(double));
    miscentered_ave_delta_sigma=(double*)malloc(Nbins*sizeof(double));
  }

  //Acquire the mass, concentration, miscentering, miscentered fraction and cosmology
  double Mass,concentration,omega_m,H0,Rmis,fmis;
  status |= c_datablock_get_double(block,cosmo,"log10_cluster_mass",&Mass);
  status |= c_datablock_get_double(block,cosmo,"omega_m",&omega_m);
  status |= c_datablock_get_double(block,cosmo,"hubble",&H0);
  status |= c_datablock_get_double(block,cosmo,"sigma_miscentering",&Rmis);
  status |= c_datablock_get_double(block,cosmo,"miscentered_fraction",&fmis);
  status |= c_datablock_get_double(block,"mc_relation","concentration",&concentration);
  Mass = pow(10,Mass);

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=200;
  params->Rmis=Rmis;
  params->fmis=fmis;
  params->timing=1; //1 is true
  params->miscentering=1;//1 is true
  params->averaging=1; //1 is true
  params->Nbins=Nbins;
  params->R_bin_min=bin_min;
  params->R_bin_max=bin_max;

  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  
  //Calculate everything one redshift bin at a time
  for( i = 0; i < NZ; i ++){
    outputs->R=R[i];
    outputs->xi_1halo=xi_nfw[i];
    outputs->xi_mm=xi_mm[i];
    outputs->xi_2halo=xi_2halo[i];
    outputs->xi_hm=xi_hm[i];
    outputs->sigma_r=sigma_r[i];
    outputs->delta_sigma=delta_sigma[i];
    outputs->bias=&bias[i];
    outputs->nu=&nu[i];
    outputs->miscentered_sigma_r=miscentered_sigma_r[i];
    outputs->miscentered_delta_sigma=miscentered_delta_sigma[i];
    outputs->Rbins=Rbins[i];
    outputs->ave_delta_sigma=ave_delta_sigma[i];
    outputs->miscentered_ave_delta_sigma=miscentered_ave_delta_sigma[i];

    status|=interface(k_lin,P_lin[i],N_lin,k_nl,P_nl[i],N_nl,
	    NR,Rmin,Rmax,*cosmo,params,outputs);
  }

  //Compute the full model
  double model_out[NZ][NR];
  double ave_model_out[NZ][Nbins];
  for(i = 0 i < NZ; i++){
    for(j = 0; j < NR; j++){
      model_out[i][j]=(1.-fmis)*delta_sigma[i][j]+fmis*miscentered_delta_sigma[i][j];
    }
    for(j = 0; j < Nbins; j++){
      ave_model_out[i][j]=(1.-fmis)*ave_delta_sigma[i][j]+fmis*ave_miscentered_delta_sigma[i][j];
    }
  }

  //Write everything to the block
  double R_out[NZ][NR];
  double xi_mm_out[NZ][NR];
  double xi_nfw_out[NZ][NR];
  double xi_2halo_out[NZ][NR];
  double xi_hm_out[NZ][NR];
  double sigma_r_out[NZ][NR];
  double delta_sigma_out[NZ][NR];
  double miscentered_sigma_r_out[NZ][NR];
  double miscentered_delta_sigma_out[NZ][NR];
  double Rbins_out[NZ][Nbins];
  double ave_delta_sigma_out[NZ][Nbins];
  double miscentered_ave_delta_sigma_out[NZ][Nbins];
  for(i = 0 i < NZ; i++){
    for(j = 0; j < NR; j++){
      R_out[i][j] = R[i][j];
      xi_mm_out[i][j] = xi_mm[i][j];
      xi_nfw_out[i][j] = xi_nfw[i][j];
      xi_2halo_out[i][j] = xi_2halo[i][j];
      xi_hm_out[i][j] = xi_hm[i][j];
      sigma_r_out[i][j] = sigma_r[i][j];
      delta_sigma_out[i][j] = delta_sigma[i][j];
      miscentered_sigma_r_out[i][j] = miscentered_sigma_r[i][j];
      miscentered_delta_sigma_out[i][j] = miscentered_delta_sigma[i][j];
    }
    for(j = 0; j < Nbins; j++){
      Rbins_out[i][j] = Rbins[i][j];
      ave_delta_sigma_out[i][j] = ave_delta_sigma[i][j];
      miscentered_ave_delta_sigma_out[i][j] = miscentered_ave_delta_sigma[i][j];
    }
  }
  status |= c_datablock_put_double_array(block,"R","R",(double*)R_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"xi_mm","xi_mm",(double*)xi_mm_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"xi_nfw","xi_nfw",(double*)xi_nfw_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"xi_2halo","xi_2halo",(double*)xi_2halo_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"xi_hm","xi_hm",(double*)xi_hm_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"sigma_r","sigma_r",(double*)sigma_r_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"delta_sigma","delta_sigma",(double*)delta_sigma_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"miscentered_sigma_r","miscentered_sigma_r",(double*)miscentered_sigma_r_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"miscentered_delta_sigma","miscentered_delta_sigma",(double*)miscentered_delta_sigma_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"delta_sigma_full","delta_sigma_full",(double*)model_out,ndim_out,extents);
  status |= c_datablock_put_double_array(block,"Rbins","Rbins",(double*)Rbins_out,ndim_out,extents_ave);
  status |= c_datablock_put_double_array(block,"ave_delta_sigma","ave_delta_sigma",(double*)ave_delta_sigma_out,ndim_out,extents_ave);
  status |= c_datablock_put_double_array(block,"miscentered_ave_delta_sigma","miscentered_ave_delta_sigma",(double*)miscentered_ave_delta_sigma_out,ndim_out,extents_ave);
  status |= c_datablock_put_double_array(block,"ave_delta_sigma_full","ave_delta_sigma_full",(double*)ave_model_out,ndim_out,extents_ave);
  status |= c_datablock_put_double_array(block,"bias","bias",bias,1,extents_bias);
  status |= c_datablock_put_double_array(block,"nu","nu",nu,1,extents_bias);

  //Free everything
  for(i = 0; i < NZ; i++){
    free(R[i]),free(xi_mm[i]),free(xi_nfw[i]),free(xi_2halo[i]),free(xi_hm[i]);
    free(sigma_r[i]),free(delta_sigma[i]),free(ave_delta_sigma[i]);
    free(miscentered_sigma_r[i]),free(miscentered_delta_sigma[i]);
    free(miscentered_ave_delta_sigma[i]);
  }
  free(R),free(xi_mm),free(xi_nfw),free(xi_2halo),free(xi_hm);
  free(sigma_r),free(delta_sigma),free(ave_delta_sigma);
  free(miscentered_sigma_r),free(miscentered_delta_sigma);
  free(miscentered_ave_delta_sigma);

  return 0
}

int cleanup(void*config_in){
  weak_lensing_config*config = (weak_lensing_config*)config_in;
  free(config);
  return 0;
}
