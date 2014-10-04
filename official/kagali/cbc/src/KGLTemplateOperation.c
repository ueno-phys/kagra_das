/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes functions relevant for inspiral waveform generation.
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-06-26 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */

#include <kagali/KGLStdlib.h>
#include <kagali/nrutil.h>
#include <kagali/KGLTemplateWaveform.h>
#include <kagali/KGLParameters.h>
#include <complex.h>


void KGLInnerProdFreq( //begin{proto}
    double complex *a1,
    double complex *a2,
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double complex *out /** only real part should be used*/
    ) //end{proto}
{
  // take the inner product of (a,b) in freq. domain
  int kFreq,kmin,kmax;
  double complex ctmp,csum;
  double df,maxf;
  size_t nfreq = (n/2+1);
  
  df=fs/n;
  kmin=fmin/df;
  
  //  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  //  maxf=MINd(fmax,MINd(fisco,fs/2.));  
  maxf=MINd(fmax,fs/2.);
  kmax=maxf/df;
  
  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }
  
  csum=0;
  for(kFreq=0;kFreq<nfreq;kFreq++){
    if(kFreq>=kmin && kFreq<=kmax && Sn[kFreq]>0){
      ctmp = a1[kFreq]*conj(a2[kFreq]);
      ctmp = ctmp/Sn[kFreq];
      csum += ctmp*df;
    }
  }
  *out = csum;
  
  return;
}


void KGLInnerProdTime( //begin{proto}
   double *a,
   double *b,
   int n,
   double *sum
   ) //end{proto}
{
  // take the inner product of (a,b) in time domain
  
  int k;
  double tmp=0;
  
  for(k=0;k<n;k++){
    tmp += a[k]*b[k];
  }
  *sum=tmp;
  
  return;
}


void KGLWhitening( //begin{proto}
    double complex *a,
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double complex *aw
    ) //end{proto}
{
  // whiten any FD vector using a given Sn  
  int kFreq,kmin,kmax;
  float df,maxf;
  size_t nfreq = (n/2+1);
  
  df=fs/n;
  kmin=fmin/df;

  //  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  //  maxf=MINd(fmax,MINd(fisco,fs/2.));  
  maxf=MINd(fmax,fs/2.);  
  kmax=maxf/df;
  
  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }
  
  for(kFreq=0;kFreq<nfreq;kFreq++){
    if(kFreq<kmin || kFreq>kmax){
      aw[kFreq]=0;
    } else {
      if(Sn[kFreq]<=0 || a[kFreq]==0){
        aw[kFreq]=0;
      }else{
        aw[kFreq]=a[kFreq]*pow(Sn[kFreq],-0.5);
      }
    }
  }
  
  return;
}
