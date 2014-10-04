/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *
 *  This program gives noise power spectrum density.
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-06-26 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */

#include <kagali/KGLStdlib.h>
#include <kagali/KGLNoisePSD.h>
#define NOISEPATH "/Users/kagra/official/kagali/detector/noisedata/"

// The following noise PSD are available at
// Damour, Iyer, Sathy, PRD63, 044023 (2001) 

void KGLNoiseGEO( //begin{proto}
    int    n,
    double *Sn,
    double fs
    ) //end{proto}
{
/* GEO */
  int i;
  double freq,df,f0,x;
  size_t nfreq = (n/2+1);
  
  f0=150.;
  df=fs/n;
  
  Sn[0]=0;
  for(i=1;i<nfreq;i++){
    freq=df*i;
    x=freq/f0;
    Sn[i]=pow(3.4*x,-30)+34.*pow(x,-1)
      +20*(1-pow(x,2)+0.5*pow(x,4)/(1+0.5*pow(x,2)));
    Sn[i]=Sn[i]*pow(10,-46);
    Sn[i]=Sn[i]/2.; // 2-sided -> 1-sided
  }
  
  return;

}


void KGLNoiseLIGOI( //begin{proto}
    int    n,
    double *Sn,
    double fs
    ) //end{proto}
{
/* initial LIGO */
  int i;
  double freq,df,f0,x;
  size_t nfreq = (n/2+1);
  
  f0=150.;
  df=fs/n;
  
  Sn[0]=0;
  for(i=1;i<nfreq;i++){
    freq=df*i;
    x=freq/f0;
    Sn[i]=9.0*(pow(4.49*x,-56.)+0.16*pow(x,-4.52)+0.52+0.32*pow(x,2.));
    Sn[i]=Sn[i]*pow(10,-46);
    Sn[i]=Sn[i]/2.; // 2-sided -> 1-sided
  }
  
  return;

}


// The following noise PSD are available at
// Cokelaer, PRD76, 102004 (2007) 

void KGLNoiseLIGOA( //begin{proto}
    int    n,
    double *Sn,
    double fs
    ) //end{proto}
{
/* advanced LIGO */
  int i;
  double freq,df,f0,x,x2;
  size_t nfreq = (n/2+1);
  
  f0=215.;
  df=fs/n;
  
  Sn[0]=0;
  for(i=1;i<nfreq;i++){
    freq=df*i;
    x=freq/f0;
    x2=x*x;
    Sn[i]=pow(x,-4.14)-5./x2+111.*(1. - x2 + 0.5*x2*x2)/(1. + 0.5*x);
    Sn[i]=Sn[i]*pow(10,-49);
    Sn[i]=Sn[i]/2.; // 2-sided -> 1-sided
  }
  
  return;
  
}


void KGLNoiseVIRGO( //begin{proto}
    int    n,
    double *Sn,
    double fs
    ) //end{proto}
{
/* VIRGO */
  int i;
  double freq,df,f0,x;
  size_t nfreq = (n/2+1);  
  f0=500.;
  df=fs/n;
  
  Sn[0]=0;
  for(i=1;i<nfreq;i++){
    freq=df*i;
    x=freq/f0;
    Sn[i]=pow(7.87*x,-4.8)+6./17.*pow(x,-1)+(1.+pow(x,2));
    Sn[i]=Sn[i]*10.2*pow(10,-46);
    Sn[i]=Sn[i]/2.; // 2-sided -> 1-sided
  }

  return;
  
}


void KGLReadNoiseSpectrum( //begin{proto}
    int    iDet,
    int    n,
    double *Sn,
    double fs
    ) //end{proto}
{
  double freqin[10000],Snin[10000]; // only used for KAGRA noise PSD
  size_t nfreq = (n/2+1);
  char filename[256];
  FILE *fp;
  
  if(iDet==1||iDet==2){
    if(iDet==1){
      sprintf(filename,"%s%s",NOISEPATH,"BW2009_VRSED.dat");
    }else if(iDet==2){
      sprintf(filename,"%s%s",NOISEPATH,"BW2009_VRSED.dat");
    }
    fp=fopen(filename,"r");
    double tmp1,tmp2;
    int num=0;
    while(1){
      if(!feof(fp)){
	num ++;
	fscanf(fp,"%lg %lg",&tmp1,&tmp2);
	freqin[num]=tmp1;
	Snin[num]=pow(tmp2,2);
      }else break;
    }
    fclose(fp);
    if(num>9999){
      fprintf(stderr,"Input data exceed the array size.\n");
      exit(EXIT_FAILURE);
    }
    
    int i,j;
    double freq,df;
    df=fs/n;
    
    // linear interpolation of the noise PSD
    Sn[0]=0;
    for(i=1;i<nfreq;i++){
      freq=df*i;
      for(j=1;j<num;j++){
	if(freq>freqin[j] && freq<=freqin[j+1]){
	  Sn[i]=Snin[j]+(Snin[j+1]-Snin[j])
	    *(freq-freqin[j])/(freqin[j+1]-freqin[j]);
	}
      }
    }

  }else if(iDet==3){
    KGLNoiseLIGOI(n,Sn,fs);
  }else if(iDet==4){
    KGLNoiseLIGOA(n,Sn,fs);
  } else {
    fprintf(stderr,"No noise curve is selected.\n");
    exit(EXIT_FAILURE);
  }

  return;
  
}
