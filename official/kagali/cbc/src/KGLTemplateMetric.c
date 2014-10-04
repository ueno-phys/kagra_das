/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program includes the functions necessary to calculate
 *  the mertric of a 2D mass parameter space.
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
#include <kagali/KGLTemplateMetric.h>
#include <kagali/KGLParameters.h>


void KGLNoiseMoment( //begin{proto}
   int    iq,
   int    ir,
   double *Sn,
   int    n,
   double fmin,
   double fmax,
   double fs,
   double *val) //end{proto}
{
  int kFreq,kmin,kmax;
  double vtemp;
  double f,df;
  double x;
  
  df=fs/n;
  kmin=fmin/df;
  kmax=fmax/df+1;
  
  vtemp=0;
  for(kFreq=0;kFreq<n/2;kFreq++){
    f=df*kFreq;
    x=f/fmin;
    if(kFreq>=kmin && kFreq<=kmax && Sn[kFreq]>0){
      vtemp += df/fmin*pow(x,-iq/3.)*pow(log(x),ir)/Sn[kFreq];
    }
  }
  
  *val = vtemp;
  
  return ;
}


void KGLNoiseMoment7( //begin{proto}
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double **vr
    ) //end{proto}
{
  int i,j;
  double val, val7;
  
  KGLNoiseMoment(7,0,Sn,n,fmin,fmax,fs,&val7);
  for(i=0;i<=17;i++){
    for(j=0;j<=2;j++){
      KGLNoiseMoment(i,j,Sn,n,fmin,fmax,fs,&val);
      vr[i][j]=val/val7;
    }
  }
  
  return;
}


void KGLPhaseGradExpansionCoeff( //begin{proto}
    double tau0,
    double tau3,
    double fmin,
    double **psi
    ) //end{proto}
{ 
  // cf. Babak et al., CQG 23 (2006) 5477; Eq. (3.24)
  // Extended to 3PN order
  int i,j;
  
  for(i=0;i<=8;i++){
    psi[1][i] = 0;
    psi[2][i] = 0;
  }  
  
  // Eq. (3.18)
  double PI=KGL_PI;
  double theta1, theta2;
  theta1 = 2.*PI*fmin*tau0;
  theta2 = 2.*PI*fmin*tau3;
  
  if(theta1<=0||theta2<=0){
    fprintf(stderr,"fmin,tau0,tau3 must be positive!");
    exit(1);
  }
  
  double lambda,theta,gamma;
  lambda=-1987/3080.;
  theta=-11831/9240.;
  gamma=0.5772156649015329;
  
  double **a,**b,**c;
  a=dmatrix(0,8,1,20);
  b=dmatrix(0,8,1,20);
  c=dmatrix(0,8,1,20);
  // initialization    
  for(i=0;i<=8;i++) {
    for(j=1;j<=20;j++) {
      a[i][j]=0;
      b[i][j]=0;
      c[i][j]=0;
    }
  }
  
  a[0][1]=3./5.;
  a[2][1]=11.*PI/12.;
  a[2][2]=743/2016.*pow(25/2./pow(PI,2),1/3.);
  a[3][1]=-3/2.;
  a[4][1]=617/384.*pow(PI,2);
  a[4][2]=5429/5376.*pow(25*PI/2.,1/3.);
  a[4][3]=15293365/10838016.*pow(5/4./pow(PI,4),1/3.);
  a[5][1]=-65.*PI/384.;
  a[5][2]=38645./64512.*pow(25/2./pow(PI,2),1/3.);
  a[5][3]=325./384.*PI*log(2);
  a[5][4]=-193225./64512.*pow(25/2./pow(PI,2),1/3.)*log(2);
  a[5][5]=-65.*PI/384.*log(5);
  a[5][6]=38645./64512.*pow(25/2./pow(PI,2),1/3.)*log(5);
  //  a[5][7]=-65.*PI/384.*log(f); can be discarded since b[5][7] and c[5][7] are zero. 
  //  a[5][8]=38645./64512.*pow(25/2./pow(PI,2),1/3.)*log(f); -> a[7][1]
  a[5][9]=65.*PI/384.*log(fmin);
  //  a[5][10]=-38645./64512.*pow(25/2./pow(PI,2),1/3.)*log(fmin); -> a[7][1]
  a[5][11]=65.*PI/384.*log(PI);
  a[5][12]=-38645./64512.*pow(25/2./pow(PI,2),1/3.)*log(PI);
  a[5][13]=65.*PI/384.;
  a[5][14]=-38645./64512.*pow(25/2./pow(PI,2),1/3.);
  a[5][15]=-65.*PI/384.;
  a[5][16]=38645./64512.*pow(25/2./pow(PI,2),1/3.);
  a[6][1]=-25565./27648.*pow(PI,3);
  a[6][2]=15211./73728.*pow(25/2.*pow(PI,4),1/3.);
  a[6][3]=-15335597827./260112384.*pow(5/4./PI,1/3.);
  a[6][4]=385./24.*lambda*pow(5/4./PI,1/3.);
  a[6][5]=2255./1024.*pow(5/4.*pow(PI,5),1/3.);
  a[6][6]=-25./8.;
  a[6][7]=11583231236531./320458457088./pow(PI,2);
  a[6][8]=-535.*gamma/112./pow(PI,2);
  a[6][9]=-55./8.*theta*pow(5/4./PI,1/3.);
  a[6][10]=-535./336./pow(PI,2)*log(2);
  a[6][11]=-535./336./pow(PI,2)*log(5);
  //  a[6][12]=-535./336./pow(PI,2)*log(f); -> a[8][1]
  //  a[6][13]=535./336./pow(PI,2)*log(fmin); -> a[8][1]
  a[6][14]=535./336./pow(PI,2)*log(PI);
  a[6][15]=535./336./pow(PI,2);
  a[6][16]=-535./336./pow(PI,2);
  //  The following terms correspond to
  //  the terms at log(x) and log(x)*x^(1/3) 
  // a[7][1] = a[5][8] + a[5][10]
  // a[8][1] = a[6][12] + a[6][13]
  a[7][1]=38645./64512.*pow(25/2./pow(PI,2),1/3.);
  a[8][1]=-535./336./pow(PI,2);
  
  b[0][1]=1.;
  b[2][1]=1./theta2;
  b[2][2]=1./3.*pow(theta2/theta1,2/3.);
  b[4][1]=1/pow(theta2,2);
  b[4][2]=1./3./pow(pow(theta1,2)*theta2,1/3.);
  b[4][3]=-1./3.*pow(theta2/theta1,4/3.);
  b[5][2]=-2./3.*pow(theta2/theta1,5/3.);
  b[5][4]=b[5][2];
  b[5][6]=b[5][2];
  //  b[5][8]=b[5][2]; corresponding a[][] is log(x) -> b[7][1]
  //  b[5][10]=b[5][2]; corresponding a[][] is log(x) -> b[7][1]
  b[5][12]=b[5][2];
  b[5][13]=1./theta1;
  b[5][14]=pow(theta2/theta1,5/3.)*(1.-2./3.*log(theta1));
  b[5][16]=-2./3.*log(theta2)*pow(theta2/theta1,5/3.);
  b[6][1]=1./pow(theta2,3);
  b[6][2]=1./3./pow(theta1*pow(theta2,2),2/3.);
  b[6][3]=-1./3.*pow(theta2/pow(theta1,4),1/3.);
  b[6][4]=b[6][3];
  b[6][5]=b[6][3];
  b[6][6]=-pow(theta2/theta1,2);
  b[6][7]=b[6][6];
  b[6][8]=b[6][6];
  b[6][9]=b[6][3];
  b[6][10]=b[6][6];
  b[6][11]=b[6][6];
  //  b[6][12]=b[6][6]; corresponding a[][] is log(x)*x^(1/3) -> b[8][1]
  //  b[6][13]=b[6][6]; corresponding a[][] is log(x)*x^(1/3) -> b[8][1]
  b[6][14]=b[6][6];
  b[6][15]=pow(theta2/theta1,2)*(1.-log(theta1));
  b[6][16]=-pow(theta2/theta1,2)*log(theta2);
  b[7][1]=b[5][2];
  b[8][1]=b[6][6];
  
  c[2][1]=-theta1/pow(theta2,2);
  c[2][2]=2./3.*pow(theta1/theta2,1/3.);
  c[3][1]=1.;
  c[4][1]=-2.*theta1/pow(theta2,3);
  c[4][2]=-1./3.*pow(theta1/pow(theta2,4.),1/3.);
  c[4][3]=4./3.*pow(theta2/theta1,1/3.);
  c[5][2]=5./3.*pow(theta2/theta1,2/3.);
  c[5][4]=c[5][2];
  c[5][6]=c[5][2];
  c[5][8]=c[5][2];
  c[5][10]=c[5][2];
  c[5][12]=c[5][2];
  c[5][14]=5./3.*pow(theta2/theta1,2/3.)*log(theta1);
  c[5][15]=1./theta2;
  c[5][16]=pow(theta2/theta1,2/3.)*(1.+5./3.*log(theta2));
  c[6][1]=-3.*theta1/pow(theta2,4);
  c[6][2]=-4./3.*pow(theta1/pow(theta2,7),1/3.);
  c[6][3]=1./3./pow(theta1*pow(theta2,2),1/3.);
  c[6][4]=c[6][3];
  c[6][5]=c[6][3];
  c[6][6]=2.*theta2/theta1;
  c[6][7]=c[6][6];
  c[6][8]=c[6][6];
  c[6][9]=c[6][3];
  c[6][10]=c[6][6];
  c[6][11]=c[6][6];
  c[6][12]=c[6][6];
  c[6][13]=c[6][6];
  c[6][14]=c[6][6];
  c[6][15]=2.*theta2/theta1*log(theta1);
  c[6][16]=theta2/theta1*(1.+2*log(theta2));
  c[7][1]=c[5][2];
  c[8][1]=c[6][6];
  
  for(i=0;i<=8;i++){
    for(j=1;j<=20;j++){
      psi[1][i] += a[i][j]*b[i][j];
      psi[2][i] += a[i][j]*c[i][j];
    }
  }
  
  return;
}


void KGLTemplateMetric( //begin{proto}
    double *Sn,
    int n,
    double tau0,
    double tau3,
    double fmin,
    double fmax,
    double fs,
    double **vr,
    double **gml
    ) //end{proto}
{
  // Babak et al., CQG (2006)
  int i,k,j,m,l;
  double *J,*K,*L;
  double **psi;
  double **gamma;
  double PI=KGL_PI;
  J=dvector(0,17);
  K=dvector(0,17);
  L=dvector(0,17);
  psi=dmatrix(1,2,0,8);
  gamma=dmatrix(0,2,0,2);
  
  for(i=0;i<=17;i++){
    J[i]=vr[i][0];
    K[i]=vr[i][1];
    L[i]=vr[i][2];
  }
  
  KGLPhaseGradExpansionCoeff(tau0,tau3,fmin,psi);
  
  // Eq. (3.28)
  gamma[0][0] = 1./2.*pow(2.*PI*fmin,2)*(J[1]-pow(J[4],2));
  for(m=1;m<=2;m++){
    gamma[0][m]=0;
    for(k=0;k<=6;k++){
      gamma[0][m] += 1./2.*(2.*PI*fmin)*psi[m][k]
        *(J[9-k]-J[4]*J[12-k]);
    }
    gamma[0][m] += 1./2.*(2.*PI*fmin)
      *(psi[m][7]*(K[4]-J[4]*K[7])+psi[m][8]*(K[3]-J[4]*K[6]));
  }
  for(m=1;m<=2;m++){
    for(l=1;l<=2;l++){
      gamma[m][l]=0;
      for(k=0;k<=6;k++){
        for(j=0;j<=6;j++){
          gamma[m][l] += 1./2.*psi[m][k]*psi[l][j]
            *(J[17-k-j]-J[12-k]*J[12-j]);
        }
      }
      for(k=0;k<=6;k++){
        gamma[m][l] += 1./2.*psi[m][k]
          *(psi[l][7]*(K[12-k]-J[12-k]*K[7])
            +psi[l][8]*(K[11-k]-J[12-k]*K[6]));
      }
      for(j=0;j<=6;j++){
        gamma[m][l] += 1./2.*psi[l][j]
          *(psi[m][7]*(K[12-j]-J[12-j]*K[7])
            +psi[m][8]*(K[11-j]-J[12-j]*K[6]));
      }
      gamma[m][l] += 1./2.*psi[m][7]*psi[l][7]*(L[7]-pow(K[7],2));
      gamma[m][l] += 1./2.*psi[m][8]*psi[l][8]*(L[5]-pow(K[6],2));
      gamma[m][l] += 1./2.*(psi[m][7]*psi[l][8]+psi[m][8]*psi[l][7])
        *(L[6]-K[6]*K[7]);
    }
  }
  
  for(m=1;m<=2;m++){
    for(l=1;l<=2;l++){
      gml[m][l] = gamma[m][l]-gamma[0][m]*gamma[0][l]/gamma[0][0];
      gml[m][l] = pow(2.*PI*fmin,2)*gml[m][l]; // theta1,2 -> tau0,3
    }
  }
  
  return;
}


void KGLTemplateMetricEigen( //begin{proto}
    double **gml,
    double *theta,
    double *lambda1,
    double *lambda2
    )  //end{proto}
{
  // ds2 = gml*dtau*dtau = lambda1 dX1*dX1 + lambda2 dX2*dX2      
  double PI=KGL_PI;
  double ttmp;
  double cost,sint;
  double cost2,sint2;
  double lamtmp1,lamtmp2;
  
  ttmp=atan2(gml[1][2]+gml[2][1],gml[2][2]-gml[1][1])/2.;
  cost=cos(ttmp);
  sint=sin(ttmp);
  cost2=cost*cost;
  sint2=sint*sint;
  lamtmp1=(gml[1][1]*cost2-gml[2][2]*sint2)/(cost2-sint2);
  lamtmp2=(gml[2][2]*cost2-gml[1][1]*sint2)/(cost2-sint2);
     
  if(fabs(lamtmp1) >= fabs(lamtmp2)){
    *lambda1 = lamtmp1;
    *lambda2 = lamtmp2;
    *theta = ttmp;
  }else{
    *lambda1 = lamtmp2;
    *lambda2 = lamtmp1;
    *theta = ttmp-PI/2.;
    cost=cos(ttmp);
    sint=sin(ttmp);
    cost2=cost*cost;
    sint2=sint*sint;
  }
  
  return;
}


void KGLTemplateMetricDistance( //begin{proto}
    double **gml,
    double dtau0,
    double dtau3,
    double *ds2
    ) //end{proto}
{
  double ds2tmp;
  ds2tmp=gml[1][1]*dtau0*dtau0+2.*gml[1][2]*dtau0*dtau3+gml[2][2]*dtau3*dtau3; 
  *ds2=ds2tmp;

  return;
}
