/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes the functions called to calculate
 *  the cost function and its derivative, etc for ILS.
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
#include <kagali/KGLParameters.h>


void KGLCostLS( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double *F
    ) //end{proto}
{
  int k;
  double Cos,tmp,Ftmp;
  double PI=KGL_PI;
  Ftmp=0;
  for(k=0;k<nframe;k++){
    Cos=cos(2.*PI*(f/fs)*k+phi);
    tmp=frame[k]-A*Cos;
    Ftmp+=pow(tmp,2)/nframe;
  }
  *F=Ftmp;
  return;
}


void KGLDerivCostLS( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double *dFdA,
    double *dFdf,
    double *dFdp
    ) //end{proto}
{
  int k;
  double PI=KGL_PI;
  double Cos,Sin;
  double tmp1,tmp2,tmp3,tmp4;
  double dFdAtmp,dFdftmp,dFdptmp;
  dFdAtmp=0;
  dFdftmp=0;
  dFdptmp=0;
  for(k=0;k<nframe;k++){
    Cos=cos(2.*PI*(f/fs)*k+phi);
    Sin=sin(2.*PI*(f/fs)*k+phi);
    tmp1=frame[k]-A*Cos;
    tmp2=-2.*tmp1*Cos;
    tmp3=2.*A*tmp1*Sin*(2.*PI*(1./fs)*k);
    tmp4=2.*A*tmp1*Sin;
    dFdAtmp+=tmp2/nframe;
    dFdftmp+=tmp3/nframe;
    dFdptmp+=tmp4/nframe;       
  }
  *dFdA=dFdAtmp;
  *dFdf=dFdftmp;
  *dFdp=dFdptmp;
  return;
}


void KGLHessian3D( //begin{proto}
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double **H
    ) //end{proto}
{
  int k;
  double PI=KGL_PI;
  double Cos,Sin;
  double tmp;
  double dFdA,dFdf,dFdp;
  double dFdAdf,dFdAdp,dFdfdp;
  double dFdA2,dFdf2,dFdp2;
  dFdA=0;
  dFdf=0;
  dFdp=0;
  dFdAdf=0;
  dFdAdp=0;
  dFdfdp=0;
  dFdA2=0;
  dFdf2=0;
  dFdp2=0;
  for(k=0;k<nframe;k++){
    Cos=cos(2.*PI*(f/fs)*k+phi);
    Sin=sin(2.*PI*(f/fs)*k+phi);
    tmp=frame[k]-A*Cos;
    dFdA+=-2.*tmp*Cos/nframe;
    dFdf+=2.*A*tmp*Sin*(2.*PI*(1./fs)*k)/nframe;
    dFdp+=2.*A*tmp*Sin/nframe;
    dFdAdf+=(2.*frame[k]-4.*A*Cos)*Sin*(2.*PI*(1./fs)*k)/nframe;
    dFdAdp+=(2.*frame[k]-4.*A*Cos)*Sin/nframe;
    dFdfdp+=2.*A*(A*Sin*Sin+tmp*Cos)*(2.*PI*(1./fs)*k)/nframe;
    dFdA2+=2.*Cos*Cos/nframe;
    dFdf2+=2.*A*(A*Sin*Sin+frame[k]*Cos)*pow(2.*PI*(1./fs)*k,2)/nframe;
    dFdp2+=2.*A*(A*Sin*Sin+frame[k]*Cos)/nframe;
  }
  H[1][1]=dFdA2;
  H[1][2]=dFdAdf;
  H[1][3]=dFdAdp;
  H[2][1]=H[1][2];
  H[2][2]=dFdf2;
  H[2][3]=dFdfdp;
  H[3][1]=H[1][3];
  H[3][2]=H[2][3];
  H[3][3]=dFdp2;
  
  return;
}


void KGLInvMatrix( //begin{proto}
    int n,
    double **H
    ) //end{proto}
{
  int i,j,k;
  double buf;
  double **Horg,**Hinv;
  Horg=dmatrix(1,n,1,n); // for the final check
  Hinv=dmatrix(1,n,1,n);
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      Horg[i][j]=H[i][j];
      //      printf("%f ",H[i][j]);
      if(i==j) Hinv[i][j]=1.;
      if(i!=j) Hinv[i][j]=0;
    }
    //    printf("\n");
  }
  
  for(i=1;i<=n;i++){
    buf=1/H[i][i];
    for(j=1;j<=n;j++){
      H[i][j]*=buf;
      Hinv[i][j]*=buf;
    }
    for(j=1;j<=n;j++){
      if(i!=j){
	buf=H[j][i];
	for(k=1;k<=n;k++){
	  H[j][k]-=H[i][k]*buf;
	  Hinv[j][k]-=Hinv[i][k]*buf;
	}
      }
    }
  }
  
  // check the inverse matrix is generated correctly.
  /*
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      buf=0;
      for(k=1;k<=n;k++){
	//	buf+=H[i][k]*Hinv[k][j];
	buf+=Horg[i][k]*Hinv[k][j];
      }
      printf("%d %d %f\n",i,j,buf);
    }
  }
  */
  
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      //      printf("%f ",Hinv[i][j]);
      H[i][j]=Hinv[i][j];
    }
    //    printf("\n");
  }
  
  return;
}

