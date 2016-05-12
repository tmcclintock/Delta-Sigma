/*
  
 */

#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_histogram.h>

int main(int argc, char **argv){
  int i;
  /*iterator*/
  FILE *fp;
  /*file pointer*/
  char fname[200];
  /*file name*/
  int n_samples = 10000; /*number of samples from pdf*/
  int n_bins = 100;
  /*number of bins in the histogram*/
  double sigma = 1.0; /*stddev of gaussian*/
  double x_min = -5.0; /*min of histogram range*/
  double x_max = 5.0; /*max of histogram range*/
  double sample;
  /*random number pulled from the pdf*/
  double norm;
  /*normalization*/
  /*GSL histogram structure*/
  gsl_histogram *h = gsl_histogram_alloc(n_bins);
  /*GSL random number generator type*/
  const gsl_rng_type *T;
  /*GSL random number generator*/
  gsl_rng *r;
  /*initialize rng*/
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*set the RNG generator to the default*/
  /*allocate a GSL RNG*/
  
  /*initialize the histogram to range [x_min,x_max)*/
  gsl_histogram_set_ranges_uniform(h, x_min, x_max);
  /*initialize the histogram to range [x_min,x_max)*/
  gsl_histogram_set_ranges_uniform(h, x_min, x_max);
  /*generate gaussian-distributed random numbers*/
  /*and histogram them*/
  for(i=0;i<n_samples;i++)
    {
      /*get gaussian rn*/
      sample = gsl_ran_gaussian(r, sigma);
      /*put it in the histogram*/
      gsl_histogram_increment(h, sample);
    }
  /*normalize histogram*/
  norm = gsl_histogram_sum(h) * (x_max-x_min)/((double) n_bins);
  gsl_histogram_scale(h, 1./norm);
  /*open a file to write the histogram*/
  sprintf(fname,"gaussian_histogram.txt");
  fp = fopen(fname,"w");
  fprintf(fp,"%d\n",n_bins);
  /*output histogram to a file*/
  gsl_histogram_fprintf(fp, h, "%f","%f");
  /*close the file*/
  fclose(fp);
  /*free memory*/
  gsl_rng_free(r);
  gsl_histogram_free(h);
  /*we're done*/
  return 0;
}
