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


void KGLGenerateFreqPow( //begin{proto}
    int n,
    double fs,
    double power,
    double *freqpow
    ) //end{proto}
{
  int kFreq;
  double df,f;
  size_t nfreq = (n/2+1);

  df=fs/n;
  
  freqpow[0]=0;
  for(kFreq=1;kFreq<nfreq;kFreq++){
    f=kFreq*df;
    freqpow[kFreq]=pow(f,power);
  }

  return;
}


void KGLGenerateTimePow( //begin{proto}
    int n,
    double fs,
    double power,
    double *timepow
    )  //end{proto}
{
  int kTime;
  double dt,t;
  
  dt=1./fs;
  timepow[0]=0;
  for(kTime=1;kTime<n;kTime++){
    t=kTime*dt;
    timepow[kTime]=pow(t,power);
  }
  
  return;
}


void KGLInspiralPhaseTD( //begin{proto}
    int iPN,
    double M,
    double eta,
    int n,
    double tmin,
    double tmax,
    double fs,
    double *phase
    ) //end{proto}
{
  // Ref. Luc Blanchet
  // http://www.livingreviews.org/lrr-2014-2
  // Eq. (317)
  // truncated at 3PN
  
  int kTime,kmin,kmax;
  double t,dt,theta;
  double gamma=0.5772156649015329;
  double theta0 = 1.; //  theta0 = ((n-1)/fs)*eta/(5.*TTSOLAR*M);
  // theta0 is a constant of integration that can be fixed by the initial conditions 
  // when the wave frequency enters the detector.
  
  dt=1/fs;
  kmin=(int)(tmin/dt);
  kmax=(int)(tmax/dt);
  
  for(kTime=kmin;kTime<kmax;kTime++){
    t=kTime*dt;
    if(t==0){
      phase[kTime]=0;
    }else{
      theta = t*eta/(5.*TTSOLAR*M);  
      if(iPN>=0){
	phase[kTime]=-1./eta*pow(theta,5./8.);
      }
      if(iPN>=2){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(3715./8064.+55./96.*eta)*pow(theta,-1./4.);
      }
      if(iPN>=3){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(-3./4.*KGL_PI)*pow(theta,-3./8.);
      }
      if(iPN>=4){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(9275495./14450688.+284875./258048.*eta+1855./2048.*pow(eta,2))*pow(theta,-1./2.);
      }
      if(iPN>=5){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(-38645./172032.+65./2048.*eta)*KGL_PI*log(theta/theta0)*pow(theta,-5./8.);
      }
      if(iPN>=6){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(831032450749357./57682522275840.-53./40.*pow(KGL_PI,2.)-107./56.*gamma
	     +107./448.*log(theta/256.)
	     +(-126510089885./4161798144.+2255./2048.*pow(KGL_PI,2.))*eta
	     +154565./1835008.*pow(eta,2.)-1179625./1769472.*pow(eta,3.))*pow(theta,-3./4.);
      }
      if(iPN>=7){
	phase[kTime]=phase[kTime]-1./eta*pow(theta,5./8.)
	  *(188516689/173408256+488825/516096*eta
	    -141769/516096*pow(eta,2))*KGL_PI*pow(theta,-7./8.);
      }	    
    }
    if(phase[kTime]!=phase[kTime]){
      printf("KGLInspiralPhaseTD: kTime=%d iPN=%d theta=%e phase=%e\n",kTime,iPN,theta,phase[kTime]);
      exit(EXIT_FAILURE);
    }
    phase[kTime]=2.*phase[kTime];
  }
  
  return;
}


void KGLInspiralPhaseFD( //begin{proto}
    int    iPN,
    double M,
    double eta,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double *phase
    ) //end{proto}
{
  int kFreq,kmin,kmax;
  double fisco,maxf;
  double f,df,x;
  double lambda=-1987/3080.;
  double theta=-11831/9240.;
  double gamma=0.5772156649015329;
  size_t nfreq = (n/2+1);
  
  df=fs/n;
  kmin=(int)(fmin/df);
  
  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  maxf=MINd(fmax,MINd(fisco,fs/2.));
  kmax=maxf/df;
  
  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }
  
  for(kFreq=kmin;kFreq<kmax;kFreq++){
    f=kFreq*df;
    x=pow(KGL_PI*TTSOLAR*M*f,1./3.);
    
    if(iPN>=0){
      phase[kFreq]=2*KGL_PI*f*tc+phic+3./128./eta/pow(x,5.);
    }
    if(iPN>=2){
      phase[kFreq]=phase[kFreq]+3./128./eta/pow(x,5.)
	*(3715./756.+55.*eta/9.)*pow(x,2.);
    }
    if(iPN>=3){     
      phase[kFreq]=phase[kFreq]+3./128./eta/pow(x,5.)
	*(-16*KGL_PI)*pow(x,3.);
    }
    if(iPN>=4){
      phase[kFreq]=phase[kFreq]+3./128./eta/pow(x,5.)
	*(15293365/508032.+27145*eta/504.+3085*eta*eta/72.)*pow(x,4.);
    }
    if(iPN>=5){
      phase[kFreq]=phase[kFreq]+3./128./eta/pow(x,5.)
	*(KGL_PI*(38645./756.-65./9.*eta)*(1.+3.*log(x)))*pow(x,5.);
    }
    if(iPN>=6){
      phase[kFreq]=phase[kFreq]+3./128./eta/pow(x,5.)
	*(11583231236531./4694215680.-15335597827./3048192.*eta
	    +2255/12.*eta*KGL_PI*KGL_PI+12320/9.*eta*lambda+76055/1728.*eta*eta
	    -127825/1296.*eta*eta*eta-640/3.*KGL_PI*KGL_PI-6848/21*gamma
	    -13696./21.*log(2.)-6848./21.*log(x)-1760./3.*eta*theta)*pow(x,6.);
    }
    if(phase[kFreq]!=phase[kFreq]){
      printf("KGLInspiralPhaseFD: kFreq=%d iPN=%d x=%e phase=%e\n",kFreq,iPN,x,phase[kFreq]);
      exit(EXIT_FAILURE);
    }
  }
  
  return;
}


void KGLGenerateInspiralTemplateFD( //begin{proto}
    int    iPN,
    int    npoint,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double M,
    double eta,
    double *freqpow,
    double complex *a1,
    double complex *a2
    ) //end{proto}
{
  int kFreq,kmax,kmin;
  double df,*phase,fisco,maxf;
  size_t nfreq = (npoint/2+1);
  
  phase=dvector(0,nfreq-1);
  KGLInspiralPhaseFD(iPN,M,eta,npoint,fmin,fmax,fs,tc,phic,phase);
  
  df=fs/npoint;
  kmin=(int)(fmin/df);
  
  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  maxf=MINd(fmax,MINd(fisco,fs/2.));  
  kmax=maxf/df;
  
  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }

  for(kFreq=0;kFreq<kmin;kFreq++){
    a1[kFreq]=0; 
    a2[kFreq]=0;
  }
  
  for(kFreq=kmin;kFreq<kmax;kFreq++){
    a1[kFreq]=freqpow[kFreq]
      *(cos(phase[kFreq]) + I*sin(phase[kFreq]));
    a2[kFreq]=I*a1[kFreq];
  }
  
  for(kFreq=kmax;kFreq<nfreq;kFreq++){
    a1[kFreq]=0; 
    a2[kFreq]=0; 
  }
  
  free_dvector(phase,0,nfreq-1);
  
  return;
}


void KGLGenerateInspiralSignalTD( //begin{proto}
    int iPN,
    int npoint,
    double tmin,
    double tmax,
    double fs,
    double M,
    double eta,
    double *timepow,
    double *a1,
    double *a2
    ) //end{proto}
{
  int kTime,kmax,kmin;
  double dt,*phase;
  
  // initialize phase
  phase=dvector(0,npoint-1);
  KGLInspiralPhaseTD(iPN,M,eta,npoint,tmin,tmax,fs,phase);
  
  dt=1./fs;
  kmin=(int)(tmin/dt);
  kmax=(int)(tmax/dt);
  
  for(kTime=0;kTime<kmin;kTime++){
    a1[kTime]=0;
    a2[kTime]=0;
  }
  for(kTime=kmin;kTime<kmax;kTime++){
    a1[kTime]=timepow[kTime]*cos(phase[kTime]);
    a2[kTime]=timepow[kTime]*sin(phase[kTime]);
  }
  for(kTime=kmax;kTime<npoint;kTime++){ 
    a1[kTime]=0;
    a2[kTime]=0;
  }
  
  free_dvector(phase,0,npoint-1);
  
  return;
}


/* freqpow[0..n-1]= f^(-7/6) */
void KGLGenerateInspiralSignalFD( //begin{proto}
    int iPN,
    int npoint,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double M,
    double eta,
    double *freqpow,
    double complex *a
    )  //end{proto}
{
  int kFreq,kmax,kmin;
  double df,*phase,fisco,maxf;
  size_t nfreq = (npoint/2+1);
  
  phase=dvector(0,nfreq-1);
  KGLInspiralPhaseFD(iPN,M,eta,npoint,fmin,fmax,fs,tc,phic,phase);
  
  df=fs/npoint;
  kmin=(int)(fmin/df);
  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  maxf=MINd(fmax,MINd(fisco,fs/2.));
  kmax=maxf/df;
  
  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }
  
  for(kFreq=0;kFreq<kmin;kFreq++){
    a[kFreq]=0;
  }
  
  for(kFreq=kmin;kFreq<kmax;kFreq++){
    a[kFreq]=freqpow[kFreq]*cos(phase[kFreq])
      +I*freqpow[kFreq]*sin(phase[kFreq]);
  }
  
  for(kFreq=kmax;kFreq<nfreq;kFreq++){ 
    a[kFreq]=0; 
  }
  
  free_dvector(phase,0,nfreq-1);
  
  return;
}


/* freqpow[0..n]= f^(-7/3) */
void KGLNormInspiral( //begin{proto}
    int npoint,
    double *Sn,
    double *freqpow,
    double fmin,
    double fmax,
    double fs,
    double M,
    double eta,
    double *norm
    ) //end{proto}
{
  int kFreq,kmin,kmax;
  double nnorm,df,fisco,maxf;
  size_t nfreq = (npoint/2+1);
  
  df=fs/npoint;
  kmin=fmin/df;

  fisco=1./(pow(6,1.5)*KGL_PI*M*TTSOLAR);
  maxf=MINd(fmax,MINd(fisco,fs/2.));
  kmax=maxf/df;

  if(kmax>nfreq){
    fprintf(stderr,"kmax=%d nfreq=%zu\n",kmax,nfreq);
    fprintf(stderr,"kmax should not exceed nfreq.\n");
    exit(EXIT_FAILURE);
  }
  
  nnorm=0.5*(freqpow[kmin]/Sn[kmin]+freqpow[kmax]/Sn[kmax]);
  for(kFreq=kmin+1;kFreq<kmax;kFreq++){
    nnorm+=freqpow[kFreq]/Sn[kFreq];
  }
  nnorm=nnorm*df;
  *norm=pow(4*nnorm,-0.5);

  return;
}


void KGLMcEtaToTau03( //begin{proto}
    double Mc,
    double eta,
    double f0,
    double *tau0,
    double *tau3
    ) //end{proto}
{
  double Mt;
  // Mukhopadhyay et al., PRD 80, 123019 (2009).
  // f0: starting frequency [Hz]
  // Mc: chirp mass
  // eta: symmetric mass ratio
  
  Mt=Mc*pow(eta,-3/5.);
  *tau0=5./(256.*KGL_PI*f0*eta)*pow(KGL_PI*Mt*f0*SOLARTIME,-5/3.);
  *tau3=1./(8.*f0*eta)*pow(KGL_PI*Mt*f0*SOLARTIME,-2/3.);

  return;
}


void KGLTau03ToMcEta( //begin{proto}
    double tau0,
    double tau3,
    double f0,
    double *Mc,
    double *eta
    ) //end{proto}
{
  double tmp1,tmp2;
  // Mukhopadhyay et al., PRD 80, 123019 (2009).
  // f0: starting frequency [Hz]
  // Mc: chirp mass
  // eta: symmetric mass ratio
  
  tmp1=pow(KGL_PI*f0,-8./5);
  tmp2=pow(KGL_PI*f0,-25./9);
  *Mc=pow(5./256/tau0,3./5)*tmp1/SOLARTIME;
  *eta=pow(KGL_PI/8./tau3,5./3)*pow(*Mc*SOLARTIME,-10./9)*tmp2;
  
  return;
}


void KGLTau03ToM12( //begin{proto}
    double tau0,
    double tau3,
    double f0,
    double *m1,
    double *m2
    ) //end{proto}
{
  double tmp1,tmp2;
  double Mc,eta,eta2;
  
  tmp1=pow(KGL_PI*f0,-8./5);
  tmp2=pow(KGL_PI*f0,-25./9);
  Mc=pow(5./256/tau0,3./5)*tmp1/SOLARTIME;
  eta=pow(KGL_PI/8./tau3,5./3)*pow(Mc*SOLARTIME,-10./9)*tmp2;
  eta2=1.-4.*eta;
  if(eta2<0 && eta2>-1e-14){
    eta2=0;
  }
  if(eta2>=0){
    *m1=(Mc/pow(eta,3./5)/2.)*(1.-pow(eta2,0.5));
    *m2=(Mc/pow(eta,3./5)/2.)*(1.+pow(eta2,0.5));
  }else{
    *m1=0;
    *m2=0;
  }

  return;
}


void KGLTime2Freq3PN( //begin{proto}
    double M,
    double eta,
    double time,
    double *freq
    ) //end{proto}
{
  // Ref. Luc Blanchet
  // http://www.livingreviews.org/lrr-2006-4
  // Eq. (233) with (232), (192) 
  // truncated at 3PN
  double x,theta;
  double gamma=0.5772156649015329;
  theta = time*eta/(5.*TTSOLAR*M);
  
  x=1./4.*pow(theta,-1./4.)
    *(1.+(743./4032.+11./48.*eta)*pow(theta,-1./4.)-1./5.*KGL_PI*pow(theta,-3./8.)
      +(19583./254016.+24401./193536.*eta+31./288.*pow(eta,2))*pow(theta,-1./2.)
      +(-11891./53760.+109./1920.*eta)*KGL_PI*pow(theta,-5./8.)
      +(-10052469856691./6008596070400.+1./6.*pow(KGL_PI,2.)+107./420.*gamma
        -107./3360.*log(theta/256.)
        +(3147553127./780337152.-451./3072.*pow(KGL_PI,2.))*eta
        -15211./442368.*pow(eta,2.)+25565./331776.*pow(eta,3.))*pow(theta,-3./4.));  
  
  if(x<0){
    *freq=0;
  }else{
    *freq=pow(x,3./2.)/(KGL_PI*TTSOLAR*M); // x=(PI*M*TTSOLAR*freq)^{2/3}
  }
  
  return;
}
