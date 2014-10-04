/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes low-pass and high-pass filters.
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-08-05 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */


#include <kagali/KGLStdlib.h>
#include <math.h>
#include <fcntl.h> 
#include <gsl/gsl_sf_bessel.h>
#include <kagali/KGLParameters.h>
#include <kagali/nrutil.h>

// low pass filter
int lpf_norm(int imode,int nkernel,double fc,double *h)
{
  //  (imode based on "Window function" of Wikipedia)
  // imode = 211 : Rectangular window
  //                    stopband attenuation -21dB (8.9%)
  //                    2.5 times faster roll-off than Blackman
  //         212 : Triangular (Bartlett) window
  //                    stopband attenuation -25dB (5.6%)
  //         231 : Hann (Hanning) window
  //                    stopband attenuation -44dB (0.63%)
  //         232 : Hamming window
  //                    stopband attenuation -53dB (0.2%)
  //         241 : Blackman window
  //                    stopband attenuation -74dB (0.02%)
  //         242 : Nuttall window
  //         243 : Blackman-Nuttall window
  //         244 : Blackman-Harris window
  //         268 : Kaiser window
  //
  // nkernel must be odd (1,3,5,...).
  // fc is normalized by fs
  //
  // transition bandwidth : delta*fc ~ 4/nkernel

  int i;
  double PI=KGL_PI;
  double hsum=0;
  double a0,a1,a2,a3;
  double alpha,beta,w;
  double arg1,tmp1,tmp2;
  
  if(fc>0.5){
    fprintf(stderr,"fc must be in the range [0,0.5].\n");
    exit(1);
  }
  
  for(i=0;i<nkernel;i++){
    if(i==(nkernel-1)/2){
      h[i]=2.*PI*fc;
    }else{
      h[i]=sin(2.*PI*fc*(i-(nkernel-1)/2))/(i-(nkernel-1)/2);
      //      printf("%d h[i]=%e\n",i,h[i]);
    }
    w=1.;
    if(imode==211){
      w=1.;
    }
    if(imode==212){
      w=1.-fabs(1.-2*i/(nkernel-1));
    }
    if(imode==231){
      w=0.5*(1.-cos(2*PI*i/(nkernel-1)));
    }
    if(imode==232){
      alpha=0.54;
      beta=1-alpha;
      w=alpha-beta*cos(2*PI*i/(nkernel-1));
    }
    if(imode==241){
      alpha=0.16;
      a0=(1-alpha)/2.;
      a1=0.5;
      a2=alpha/2.;
      w=a0-a1*cos(2*PI*i/(nkernel-1))+a2*cos(4*PI*i/(nkernel-1));
    }
    if(imode==242){
      a0=0.355768;
      a1=0.487396;
      a2=0.144232;
      a3=0.012604;
      w=a0-a1*cos(2*PI*i/(nkernel-1))+a2*cos(4*PI*i/(nkernel-1))-a3*cos(6*PI*i/(nkernel-1));
    }
    if(imode==243){
      a0=0.3635819;
      a1=0.4891775;
      a2=0.1365995;
      a3=0.0106411;
      w=a0-a1*cos(2*PI*i/(nkernel-1))+a2*cos(4*PI*i/(nkernel-1))-a3*cos(6*PI*i/(nkernel-1));
    }
    if(imode==244){
      a0=0.35875;
      a1=0.48829;
      a2=0.14128;
      a3=0.01168;
      w=a0-a1*cos(2*PI*i/(nkernel-1))+a2*cos(4*PI*i/(nkernel-1))-a3*cos(6*PI*i/(nkernel-1));
    }
    if(imode==268){
      alpha=(nkernel-1)/2;
      beta=6.;
      arg1=beta*sqrt(1.-pow((i-alpha)/alpha,2));
      tmp1=gsl_sf_bessel_I0(arg1);
      tmp2=gsl_sf_bessel_I0(beta);
      //      printf("i=%d tmp1=%f tmp2=%f\n",i,tmp1,tmp2);
      w=tmp1/tmp2; // Actually tmp2 is not necessary since h[i] is normalized later.
    }
    h[i]=h[i]*w;
    hsum += h[i];
  }
  for(i=0;i<nkernel;i++){
    h[i]=h[i]/hsum;
  }
  
  return 0;
}


// high pass filter
int hpf_norm(int imode,int nkernel,double fc,double *h)
{
  int i;
  
  if(fc>0.5){
    fprintf(stderr,"fc must be in the range [0,0.5].\n");
    exit(EXIT_FAILURE);
  }
  
  lpf_norm(imode,nkernel,fc,h);
  for(i=0;i<nkernel;i++){
    h[i]=-h[i];
    //      printf("%d h[i]=%e\n",i,h[i]);
  }
  h[(nkernel-1)/2]++;
  
  return 0;
}
