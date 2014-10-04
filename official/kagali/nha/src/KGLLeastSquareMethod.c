/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes the algorithm of 3D (A,f,phi) Newtonian's method 
 *  to obtain the minimum of a given cost function.
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
#include <kagali/KGLTemplateWaveform.h>
#include <kagali/KGLLeastSquareFunc.h>

void KGLSteepestDescent( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double mu,
    double A,
    double *f,
    double *phi,
    double *Fcost
    ) //end{proto}
{
  int i,k;
  double A0;
  double ftmp,ptmp;
  double Fcosttmp,Fcost0;
  double dFdA,dFdf,dFdp;
  double *frametmp;
  
  ftmp=*f;
  ptmp=*phi;
  
  A0=A;
  //  printf("===Steepest descent method===\n");
  frametmp=dvector(0,nframe);
  for(k=0;k<nframe;k++){
    frametmp[k]=frame[k]/A0;
    //    frametmp[k]=frame[k];
  }
  A=1;
  for(i=1;i<=nitr;i++){
    KGLCostLS(fs,nframe,frametmp,A,ftmp,ptmp,&Fcosttmp);
    KGLDerivCostLS(fs,nframe,frametmp,A,ftmp,ptmp,&dFdA,&dFdf,&dFdp);
    if(i==1) Fcost0=Fcosttmp;
    /*
    printf("i=%d A=%f f=%f phi=%f Fcost=%f\n",i,A,ftmp,ptmp,Fcost);
    printf("mu=%f dFdf=%f dFdp=%f\n",mu,dFdf,dFdp);
    */
    ftmp=ftmp-mu*dFdf;
    ptmp=ptmp-mu*dFdp;
    if(Fcosttmp<(1.-0.5*mu)*Fcost0) break;
    mu=0.5*mu;
  }
  A=A0;
  *f=ftmp;
  *phi=ptmp;
  *Fcost=Fcosttmp;
  //  printf("A=%f f=%f phi=%f Fcost=%f\n",A,ftmp,ptmp,Fcosttmp);
  
  return;
}

void KGLNewton3D( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double mu,
    double *A,
    double *f,
    double *phi
    ) //end{proto}
{
  int i;
  double Atmp,ftmp,ptmp;
  double Atmp2,ftmp2,ptmp2;
  double dF,Fcost,Fcostpre;
  double dFdA,dFdf,dFdp;
  double dA,df,dphi;
  double **H;
  H=dmatrix(1,3,1,3);
  Atmp=*A;
  ftmp=*f;
  ptmp=*phi;
  dFdA=0;
  dFdf=0;
  dFdp=0;
  KGLCostLS(fs,nframe,frame,Atmp,ftmp,ptmp,&Fcostpre);
  dF=0;
  //  printf("===Newton's method===\n");
  for(i=1;i<=nitr;i++){
    if(dF>0){
      KGLCostLS(fs,nframe,frame,Atmp,ftmp,ptmp,&Fcostpre);
    }
    //    printf("i=%d A=%f f=%f phi=%f Fcostpre=%e mu=%e\n",i,Atmp,ftmp,ptmp,Fcostpre,mu);
    KGLHessian3D(fs,nframe,frame,Atmp,ftmp,ptmp,H);
    KGLInvMatrix(3,H);
    KGLDerivCostLS(fs,nframe,frame,Atmp,ftmp,ptmp,&dFdA,&dFdf,&dFdp);
    dA=H[1][1]*dFdA+H[1][2]*dFdf+H[1][3]*dFdp;
    df=H[2][1]*dFdA+H[2][2]*dFdf+H[2][3]*dFdp;
    dphi=H[3][1]*dFdA+H[3][2]*dFdf+H[3][3]*dFdp;
    Atmp2=Atmp-mu*dA;
    ftmp2=ftmp-mu*df;
    ptmp2=ptmp-mu*dphi;
    KGLCostLS(fs,nframe,frame,Atmp2,ftmp2,ptmp2,&Fcost);
    dF=Fcostpre-Fcost;
    if(dF<0){
      /*
      printf("dF<0: i=%d Fcost=%e dF=%e\n",i,Fcost,dF);
      printf("dA=%f df=%f dphi=%f\n",dA,df,dphi);
      */
      mu=0.5*mu;
    }else{
      Atmp=Atmp2;
      ftmp=ftmp2;
      ptmp=ptmp2;
      if(Atmp<1.) Atmp=1.;
    }
  }
  *A=Atmp;
  *f=ftmp;
  *phi=ptmp;
  //  printf("A=%f f=%f phi=%f\n",Atmp,ftmp,ptmp);
  
  return;
}


void KGLAmpConv( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double nu,
    double f,
    double phi,
    double *A,
    double *Fcost
    ) //end{proto}
{
  int i;
  double Fcosttmp,Fcost0;
  double dFdA,dFdf,dFdp;
  double Atmp;
  
  Atmp=*A;
  
  //  printf("===Amplitude Convergence===\n");
  for(i=1;i<=nitr;i++){
    KGLCostLS(fs,nframe,frame,Atmp,f,phi,&Fcosttmp);
    KGLDerivCostLS(fs,nframe,frame,Atmp,f,phi,&dFdA,&dFdf,&dFdp);
    if(i==1) Fcost0=Fcosttmp;
    /*
    printf("i=%d A=%f f=%f phi=%f Fcost=%e\n",
	   i,Atmp,f,phi,Fcosttmp);
    */
    //    printf("nu=%f dFdA=%f\n",nu,dFdA);
    
    Atmp=Atmp-nu*dFdA;
    if(Fcosttmp<(1.-0.5*nu)*Fcost0) break;
    nu=0.5*nu;
  }
  *A=Atmp;
  *Fcost=Fcosttmp;
  /*
  printf("A=%f f=%f phi=%f Fcost=%e\n",Atmp,f,phi,Fcosttmp);
  printf("\n");
  */
  
  return;
}
