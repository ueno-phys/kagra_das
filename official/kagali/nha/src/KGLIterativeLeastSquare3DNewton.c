/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program iteratively least-square-fits given time-series data
 *  with a sinusoidal wave using a cost function.
 *  In other words, this program can extract any number of sinusoidal waves
 *  buried in a given time-series data in descending order of amplitude.
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-08-05 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */

#include <kagali/KGLStdlib.h>
#include <kagali/nrutil.h>
#include <kagali/RealFFT.h>
#include <kagali/KGLIterativeLeastSquare3DNewton.h>
#include <kagali/KGLLeastSquareMethod.h>
#include <kagali/KGLLeastSquareFunc.h>
#include <kagali/KGLParameters.h>

typedef double complex cdouble;

#define CALLOC(xx,nn,ss) do {				    \
    if((xx = (ss *)calloc((nn),sizeof(ss))) == NULL) {      \
      KGLPrintError("fail to alloc \"%s\"",#xx);	    \
      exit(EXIT_FAILURE);				    \
    }                                                       \
  } while(0)


void KGLIterativeLeastSquare3DNewton( //begin{proto}
   double *frame,
   double fs,
   int nframe,
   double Fcostthr,
   double a2thr,
   int nsig,
   int nitr,
   double mu0,
   double nu0,
   double *Afit,
   double *pfit,
   double *ffit
   ) //end{proto}
{
  int i,j,k;
  int kmax;
  int isig,isigpre; 
  double PI=KGL_PI;
  double A0,A,f,phi,phimin; 
  double Fcost; // cost function
  double Fcostpre,dF;
  double Fcostmin;
  double mu;
  double fmax,df_frame;
  double dA,df,dphi;
  double Atmp,ftmp,ptmp;
  double **H;
  double dFdA,dFdf,dFdp;
  double *data_vr, *data_vi;
  double **spec;
  double *power,powermax;
  double powersum;
  
  KGLStatus *status = KGLCreateStatus();
  KGLRealFFTPlan *plan = NULL;  
  size_t nframef = (nframe/2+1);
  cdouble *framef = NULL;
  CALLOC(framef,nframef,cdouble);
  
  df_frame=fs/nframe;
  data_vr=dvector(0,nframe/2);
  data_vi=dvector(0,nframe/2);
  spec=dmatrix(0,nsig,0,nframe-1);
  power=dvector(0,nframe/2);  
  for(k=0;k<nframe;k++){
    spec[0][k]=0;
  }
  H=dmatrix(1,3,1,3);
  for(isig=0;isig<nsig;isig++){
    for(k=0;k<nframe;k++){
      frame[k]=frame[k]-spec[isig][k];
    }
    
// First we find the frequency which gives the largest Amplitude (fmax) via DFT
    KGLForwardRealFFT(status,&plan,framef,frame,nframe);
    KGLAbortIfError(status);
    
    powermax=0;
    for(k=0;k<=nframe/2;k++){
      power[k]=pow(cabs(framef[k]),2);
      //	printf("%e %e\n",k*fs/nframe,power[k]);
      if(powermax<power[k]){
	powermax=power[k];
	kmax=k;
	fmax=k*fs/nframe;
      }
    }
    //  printf("arg=%f\n",acos(dataf[2*kmax]/sqrt(powermax))/PI);
    if(kmax==0){
      kmax=1;
      fmax=fs/nframe;
    }
    powersum=0;
    for(k=0;k<=nframe/2;k++){
      if(abs(kmax-k)<=2){ // freq. resolution should be taken into account
	powersum+=power[k];
      }
    }
    
    A0=2.*sqrt(powersum)/nframe;
    A=A0;
    phi=atan2(framef[2*kmax+1],framef[2*kmax]);
    f=fmax;
    //      printf("Amp=%f fmax=%f phi=%f\n",A,f,phi);      
    Fcost=0;
    Fcostpre=0;
    mu=mu0;
    dFdA=0;
    dFdf=0;
    dFdp=0;
    KGLCostLS(fs,nframe,frame,A,f,phi,&Fcostpre);
    dF=0;
    
    for(i=1;i<=nitr;i++){
      if(dF>0){
	KGLCostLS(fs,nframe,frame,A,f,phi,&Fcostpre);
      }
      //              printf("i=%d A=%f f=%f phi=%f Fcostpre=%e mu=%e\n",i,A,f,phi,Fcostpre,mu);
      KGLHessian3D(fs,nframe,frame,A,f,phi,H);
      KGLInvMatrix(3,H);
      KGLDerivCostLS(fs,nframe,frame,A,f,phi,&dFdA,&dFdf,&dFdp);
      dA=H[1][1]*dFdA+H[1][2]*dFdf+H[1][3]*dFdp;
      df=H[2][1]*dFdA+H[2][2]*dFdf+H[2][3]*dFdp;
      dphi=H[3][1]*dFdA+H[3][2]*dFdf+H[3][3]*dFdp;
      Atmp=A-mu*dA;
      ftmp=f-mu*df;
      ptmp=phi-mu*dphi;
      KGLCostLS(fs,nframe,frame,Atmp,ftmp,ptmp,&Fcost);
      dF=Fcostpre-Fcost;
      if(dF<0){
	//            printf("dF<0: i=%d Fcost=%e dF=%e\n",i,Fcost,dF);
	//            printf("dA=%f df=%f dphi=%f\n",dA,df,dphi);
	mu=0.5*mu;
      }else{
	A=Atmp;
	f=ftmp;
	phi=ptmp;       
	if(fabs(A-A0)>A0*0.1){
	  //      printf("A=%f: Amplitude is too away!!!\n",A);
	  A=A0;
	  mu=1.;
	}
	if(fabs(fmax-f)>df_frame){
	  //      printf("f=%f: Frequency is too away!!!\n",f);
	  f=fmax;
	  mu=1.;
	}
      }
      if(dF>0 && dF<Fcostthr) break;
      if(i%100==0){
	KGLCostLS(fs,nframe,frame,A,f,phi,&Fcost);
	Fcostmin=Fcost;
	phimin=phi;
	//      printf("i=%d phi=%f Fcost=%e\n",i,phi,Fcost);
	for(j=1;j<=10;j++){
	  phi=PI*j/10.;
	  KGLCostLS(fs,nframe,frame,A,f,phi,&Fcost);
	  //      printf("phi=%f Fcost=%e\n",phi,Fcost);
	  if(Fcostmin>Fcost){
	    Fcostmin=Fcost;
	    phimin=phi;
	  }
	}
	phi=phimin;
      }
    }
    if(isig<nsig){
      for(k=0;k<nframe;k++){
	spec[isig+1][k]=A*cos(2.*PI*k*(f/fs)+phi);
      }
    }
    if(A<0){
      A=-A;
      phi=phi+PI;
    }
    if(f<0){
      f=-f;
      phi=-phi;
    }
    phi=phi-2*PI*(int)(phi/2./PI);
    if(phi>PI) phi=phi-2.*PI;
    if(phi<-PI) phi=phi+2.*PI;
    if(A<a2thr) break;
    /* 
       printf("isig=%d A=%f f=%f phi=%f Fcost=%e\n",
       isig,A,f,phi,Fcost);
    */
    Afit[isig]=A;
    pfit[isig]=phi;
    ffit[isig]=f;
    if(Fcost<Fcostthr) break;
    isigpre=isig;
  }
  
  free(framef);
  free_dvector(data_vr,0,nframe/2);
  free_dvector(data_vi,0,nframe/2);
  free_dvector(power,0,nframe/2);
  free_dmatrix(spec,0,nsig,0,nframe-1);
  free_dmatrix(H,1,3,1,3);
  
  return;
}
