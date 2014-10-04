/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *
 *  This program generates a template bank based on the hexagonal placement except
 *  around the right and left edges, where a slightly different placement is applied.
 *
 *
 *  Input:
 *           iDet      (int)      assumed detector
 *                                                         1: KAGRA vRSED
 *                                                         2: KAGRA vRSEB
 *                                                         3: initial LIGO
 *                                                         4: advanced LIGO
 *
 *           fmin      (double)   minimun frequency [Hz]
 *           fmax      (double)   maximum frequency [Hz]
 *           fs        (double)   sampling rate [Hz]
 *           npoint    (int)      number of data points
 *           mm        (double)   minimal match
 *           m1min     (double)   lower threshold of m1
 *           m1max     (double)   upper threshold of m1
 *           m2min     (double)   lower threshold of m2
 *           m2max     (double)   upper threshold of m2
 *           m1ini     (double)   initial m1 value (should not be in the corner)
 *           m2ini     (double)   initial m2 value (should not be in the corner)
 *           nGridMax  (int)      maximum number of grids       (mainly used for test)
 *           nGenMax   (int)      maximum number of generations (mainly used for test)
 *
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-06-26 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp) 
 *
 */


#include <kagali/KGLStdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h> 

#include <kagali/nrutil.h>
#include <kagali/KGLTimeStamp.h>
#include <kagali/KGLTemplateWaveform.h>
#include <kagali/KGLParameters.h>
#include <kagali/KGLNoisePSD.h>
#include <kagali/KGLTemplateMetric.h>
#include <kagali/KGLHexagonReproduction.h>


int main(int argc,char *argv[])
{
  int npoint;
  int iDet;
  int nfile;
  int *iGridSt;
  int iGrid,iNewGrid,iPastGrid;
  int iVertex;
  int nGrid,nGrid2,nGridMax;
  int nGrid2tmp;
  int nNewGridMax=6;
  int nVertexMax=6;
  int *gridflag;
  int genflag;
  int *iGen,iGenLoop;
  int nGenMax;

  double ds2,mm,mmeff;
  double eta,Mc,eta2,etamin,etamax;
  double fmin,fmax,fs;
  double df;
  double tau0,tau3,dtau0,dtau3;
  double dtau0p,dtau0min,dtau0max;
  double dtau3p,dtau3min,dtau3max;
  double m1,m2,m1min,m1max,m2min,m2max;
  double Mt,Mtmin,Mtmax;
  double m1ini,m2ini;
  
  FILE *fout=NULL;
  clock_t clock0;
  
  if(argc!=18){
    fprintf(stderr,"argc=%d\n",argc);
    fprintf(stderr,"Usage: ./KGLGeneTemplateBankHexagonal "
	    "gridpoint.dat vertex.dat ellipse.dat "
	    "iDet fmin[Hz] fmax[Hz] fs[Hz] npoint "
	    "mm m1min m1max m2min m2max m1ini m2ini nGridMax nGenMax\n");
    exit(EXIT_FAILURE);
  }
  
  nfile = 3;
  iDet=atoi(argv[nfile+1]);
  fmin=atof(argv[nfile+2]);
  fmax=atof(argv[nfile+3]);
  fs=atof(argv[nfile+4]);
  npoint=atoi(argv[nfile+5]);
  mm=atof(argv[nfile+6]);
  m1min=atof(argv[nfile+7]);
  m1max=atof(argv[nfile+8]);
  m2min=atof(argv[nfile+9]);
  m2max=atof(argv[nfile+10]);
  m1ini=atof(argv[nfile+11]);
  m2ini=atof(argv[nfile+12]);
  nGridMax=atoi(argv[nfile+13]);
  nGenMax=atoi(argv[nfile+14]);
  
  df=fs/npoint;
  
  printf("You set...\n");
  printf("noise data: ");
  if(iDet==1) printf("KAGRA vRSE(D)\n");
  if(iDet==2) printf("KAGRA vRSE(B)\n");
  if(iDet==3) printf("initial LIGO\n");
  if(iDet==4) printf("advanced LIGO\n");
  printf("npoint   : %d\n",npoint);
  printf("fmin     : %6.1f\n",fmin);
  printf("fmax     : %6.1f\n",fmax);
  printf("fs       : %6.1f\n",fs);
  printf("df       : %7.5f\n",df);
  printf("mm       : %5.3f\n",mm);
  printf("m1min    : %5.2f\n",m1min);
  printf("m1max    : %5.2f\n",m1max);
  printf("m2min    : %5.2f\n",m2min);
  printf("m2max    : %5.2f\n",m2max);
  printf("m1ini    : %5.2f\n",m1ini);
  printf("m2ini    : %5.2f\n",m2ini);
  printf("nGridMax : %d\n",nGridMax);
  printf("nGenMax  : %d\n",nGenMax);

  if(m1min>m2min || m1max>m2max){
    fprintf(stderr,"m1 should be less than or equal to m2.\n");
    exit(EXIT_FAILURE);
  }
  
  clock0 = clock();
  KGLTimeStamp(clock0);

  double *Sn;
  Sn=dvector(0,npoint/2-1);
  KGLReadNoiseSpectrum(iDet,npoint,Sn,fs);

  double **vr;
  vr=dmatrix(0,17,0,2);
  KGLNoiseMoment7(Sn,npoint,fmin,fmax,fs,vr);
  
  KGLTimeStamp(clock0);
  
  double tau0min,tau3min;
  double tau0max,tau3max;
  double *tau0g,*tau3g;
  double **tau0v1,**tau3v1;
  double **tau0v2,**tau3v2;
  double **gml;
  double theta;
  double lambda1,lambda2;
  double *tau0d,*tau3d;
  double *tau0vtmp,*tau3vtmp;
  iGridSt=ivector(1,nGenMax);
  tau0d=dvector(1,nNewGridMax);
  tau3d=dvector(1,nNewGridMax);
  tau0vtmp=dvector(1,nVertexMax);
  tau3vtmp=dvector(1,nVertexMax);
  iGen=ivector(1,nGridMax);
  gridflag=ivector(1,nGridMax);
  tau0g=dvector(1,nGridMax);
  tau3g=dvector(1,nGridMax);
  tau0v1=dmatrix(1,nGridMax,1,nVertexMax);
  tau3v1=dmatrix(1,nGridMax,1,nVertexMax);
  tau0v2=dmatrix(1,nGridMax,1,nVertexMax);
  tau3v2=dmatrix(1,nGridMax,1,nVertexMax);
  gml=dmatrix(1,2,1,2);
  for(iGrid=1;iGrid<=nGridMax;iGrid++){
    gridflag[iGrid]=0;
    tau0g[iGrid]=0;
    tau3g[iGrid]=0;
  }
  
  Mtmin=m1min+m2min;
  Mtmax=m1max+m2max;  

  Mt=m1min+m2min;
  eta=m1min*m2min/pow(Mt,2);
  Mc=Mt*pow(eta,3/5.);
  KGLMcEtaToTau03(Mc,eta,fmin,&tau0max,&tau3max);
  
  Mt=m1max+m2max;
  eta=m1max*m2max/pow(Mt,2);
  Mc=Mt*pow(eta,3/5.);
  KGLMcEtaToTau03(Mc,eta,fmin,&tau0min,&tau3min);
  
  m1=m1ini;
  m2=m2ini;
  Mt=m1+m2;
  eta=m1*m2/pow(Mt,2);
  Mc=Mt*pow(eta,3/5.);
  KGLMcEtaToTau03(Mc,eta,fmin,&tau0,&tau3); 
  
  iGrid=1;
  nGrid=1;
  iGen[iGrid]=1;
  iGenLoop=1;
  iGridSt[iGenLoop]=1;
  tau0g[iGrid]=tau0;
  tau3g[iGrid]=tau3;
  while(1){
    iGridSt[iGenLoop+1]=nGrid+1;
    genflag=0;
    for(iGrid=1;iGrid<=nGrid;iGrid++){
      if(iGen[iGrid]==iGenLoop){
	KGLTemplateMetric(Sn,npoint,tau0g[iGrid],tau3g[iGrid],fmin,fmax,fs,vr,gml);
	KGLTemplateMetricEigen(gml,&theta,&lambda1,&lambda2);
        if(tau0g[iGrid]>0.01*tau0max){
          mmeff=mm;
        }else{
          mmeff=mm*1.01;
        }
	
	// Reproduction
        KGLTau03ToM12(tau0g[iGrid],tau3g[iGrid],fmin,&m1,&m2);
        if(m1>m1min && m1<m1max && m2>m2min && m2<m2max){
	  KGLReprodHexagonGrids(1,mmeff,lambda1,lambda2,theta,tau0g[iGrid],tau3g[iGrid],tau0d,tau3d);
	  for(iNewGrid=1;iNewGrid<=nNewGridMax;iNewGrid++){
	    // Check if the hexagon around each new grid intersects the physical parameter region.
	    // 1. Center of the hexagon
	    KGLTau03ToMcEta(tau0d[iNewGrid],tau3d[iNewGrid],fmin,&Mc,&eta);
	    eta2=1.-4.*eta;
	    KGLTau03ToM12(tau0d[iNewGrid],tau3d[iNewGrid],fmin,&m1,&m2);
	    if(m1>m1min && m1<m1max && m2>m2min && m2<m2max){
	      gridflag[nGrid+1]=1;
	    }
	    // 2. Vertices of the hexagon
	    KGLTemplateMetric(Sn,npoint,tau0d[iNewGrid],tau3d[iNewGrid],fmin,fmax,fs,vr,gml);
	    KGLTemplateMetricEigen(gml,&theta,&lambda1,&lambda2);
	    KGLReprodHexagonVertices(1,mmeff,lambda1,lambda2,theta,
				     tau0d[iNewGrid],tau3d[iNewGrid],tau0vtmp,tau3vtmp);
            for(iVertex=1;iVertex<=nVertexMax;iVertex++){
              tau0v1[iNewGrid][iVertex]=tau0vtmp[iVertex];
              tau3v1[iNewGrid][iVertex]=tau3vtmp[iVertex];
            }
	    for(iVertex=1;iVertex<=nVertexMax;iVertex++){
	      KGLTau03ToMcEta(tau0vtmp[iVertex],tau3vtmp[iVertex],fmin,&Mc,&eta);
	      eta2=1.-4.*eta;
	      KGLTau03ToM12(tau0vtmp[iVertex],tau3vtmp[iVertex],fmin,&m1,&m2);
	      if(m1>m1min && m1<m1max && m2>m2min && m2<m2max){
		gridflag[nGrid+1]=1;
	      }
	    }
	    // 3. Maximum & minimum of tau0 & tau3 of the ellipses
	    // Upper half
	    dtau0=gml[1][2]*gml[2][1]*(1.-mmeff)/gml[1][1]
	      /(gml[1][1]*gml[2][2]-gml[1][2]*gml[2][1]);
	    dtau0=sqrt(fabs(dtau0));
	    if(gml[1][2]/(gml[1][2]*gml[2][1]-gml[1][1]*gml[2][2])<0) dtau0=-dtau0;
	    dtau3p=pow(gml[1][2]*dtau0,2)-gml[2][2]*(gml[1][1]*pow(dtau0,2)-(1.-mmeff));
	    dtau3max=(-gml[1][2]*dtau0+sqrt(dtau3p))/gml[2][2];
	    KGLTau03ToMcEta(tau0d[iNewGrid]+dtau0,tau3d[iNewGrid]+dtau3max,fmin,&Mc,&etamin);
	    // Lower half
	    dtau0=-dtau0;
	    dtau3p=pow(gml[1][2]*dtau0,2)-gml[2][2]*(gml[1][1]*pow(dtau0,2)-(1.-mmeff));
	    dtau3min=(-gml[1][2]*dtau0-sqrt(dtau3p))/gml[2][2];
	    KGLTau03ToMcEta(tau0d[iNewGrid]+dtau0,tau3d[iNewGrid]+dtau3min,fmin,&Mc,&etamax);
	    // Right half
	    dtau3=gml[1][2]*gml[2][1]*(1.-mmeff)/gml[2][2]
	      /(gml[1][1]*gml[2][2]-gml[1][2]*gml[2][1]);
	    dtau3=sqrt(fabs(dtau3));
	    if(gml[1][2]/(gml[1][2]*gml[2][1]-gml[1][1]*gml[2][2])<0) dtau3=-dtau3;
	    dtau0p=pow(gml[1][2]*dtau3,2)-gml[1][1]*(gml[2][2]*pow(dtau3,2)-(1.-mmeff));
	    dtau0max=(-gml[1][2]*dtau3+sqrt(dtau0p))/gml[1][1];
	    // Left half
	    dtau3=-dtau3;
	    dtau0p=pow(gml[1][2]*dtau3,2)-gml[1][1]*(gml[2][2]*pow(dtau3,2)-(1.-mmeff));
	    dtau0min=(-gml[1][2]*dtau3-sqrt(dtau0p))/gml[1][1];
	    if(etamin<0.25 && etamax>0.25 &&
	       tau0d[iNewGrid]+dtau0max>tau0min && tau0d[iNewGrid]-dtau0min<tau0max){
	      gridflag[nGrid+1]=1;
	    }
	    
	    // Check if the grid is produced for the first time or not.
	    if(iGenLoop<3){
	      for(iPastGrid=1;iPastGrid<=nGrid;iPastGrid++){
		KGLTemplateMetricDistance(gml,tau0g[iPastGrid]-tau0d[iNewGrid],
					  tau3g[iPastGrid]-tau3d[iNewGrid],&ds2);
		if(ds2<(1.-mmeff)/2.) gridflag[nGrid+1]=0;
	      }
	    }else{
	      for(iPastGrid=iGridSt[iGenLoop-2];iPastGrid<=nGrid;iPastGrid++){
		KGLTemplateMetricDistance(gml,tau0g[iPastGrid]-tau0d[iNewGrid],
					  tau3g[iPastGrid]-tau3d[iNewGrid],&ds2);
		if(ds2<(1.-mmeff)*0.8) gridflag[nGrid+1]=0;
	      }
	    }
	    if(gridflag[nGrid+1]==1){
	      genflag=1;
	      nGrid++;
	      iGen[nGrid]=iGenLoop+1;
	      tau0g[nGrid]=tau0d[iNewGrid];
	      tau3g[nGrid]=tau3d[iNewGrid];
	      if(nGrid % 2000 ==0){
		KGLTau03ToM12(tau0g[nGrid],tau3g[nGrid],fmin,&m1,&m2);
		printf("iGenLoop=%5d iGrid=%6d nGrid=%6d tau0=%7.2f tau3=%5.2f m1=%5.2f m2=%5.2f\n",
		       iGenLoop,iGrid,nGrid,tau0g[nGrid],tau3g[nGrid],m1,m2);
	      }
	    }
	  }
	}
	// end of reproduction
      }else if(iGen[iGrid]>iGenLoop){
	iGenLoop++;
	break;
      }
    }
    if(genflag==0) break;
    if(nGrid>=nGridMax) break;
  }
  
  printf("Then, fill the smallest and largest mass regions...\n");
// The smallest and largest mass regions are filled along eta=1/4 line.
  double A0,A3,r03;
  double tau0gmin=9999.,tau3gmin;
  double tau0gmax=0,tau3gmax;
  A0=5./256./pow(KGL_PI*fmin,8/3.);
  A3=KGL_PI/8./pow(KGL_PI*fmin,5/3.);
  
  for(iGrid=1;iGrid<=nGrid;iGrid++){
    if(tau0gmin>tau0g[iGrid]){
      tau0gmin=tau0g[iGrid];
      tau3gmin=tau3g[iGrid];
    }
    if(tau0gmax<tau0g[iGrid]){
      tau0gmax=tau0g[iGrid];
      tau3gmax=tau3g[iGrid];
    }
  }
  // First, the smallest mass region
  tau0=tau0gmax;
  tau3=4.*A3*pow(tau0/4./A0,2/5.);
  KGLTemplateMetric(Sn,npoint,tau0,tau3,fmin,fmax,fs,vr,gml);
  dtau0=(1-mm)/gml[1][1];
  dtau0=sqrt(fabs(dtau0));
  tau0=tau0-2.*dtau0;
  KGLTau03ToMcEta(tau0,tau3,fmin,&Mc,&eta);
  Mt=Mc*pow(eta,-3/5.);
  nGrid2=nGrid+1;
  tau0g[nGrid2]=tau0;
  tau3g[nGrid2]=4.*A3*pow(tau0/4./A0,2/5.);
  while(Mt>Mtmin){
    if(nGrid2 % 500 ==0){
      printf("nGrid=%6d tau0=%7.2f tau3=%5.2f\n",nGrid2,tau0,tau3);
    }
    KGLTemplateMetric(Sn,npoint,tau0,tau3,fmin,fmax,fs,vr,gml);
    
    // dtau0=r03*dtau3 on eta=1/4.
    r03=5/2.*pow(tau3,3/2.)*4*A0/pow(4*A3,5/2.);
    dtau3=(1-mm)/(gml[1][1]*pow(r03,2)+2.*gml[1][2]*r03+gml[2][2]);
    dtau3=sqrt(fabs(dtau3)); // plus sign: to increase tau0
    dtau0=r03*dtau3;
    
    tau0=tau0+1.5*dtau0;
    tau3=4.*A3*pow(tau0/4./A0,2/5.);
    KGLTau03ToMcEta(tau0,tau3,fmin,&Mc,&eta);
    Mt=Mc*pow(eta,-3/5.);
    
    nGrid2++;
    tau0g[nGrid2]=tau0;
    tau3g[nGrid2]=tau3;
  }
  
  nGrid2tmp=nGrid2;
  for(iGrid=nGrid+1;iGrid<=nGrid2tmp;iGrid++){
    tau0=tau0g[iGrid];
    tau3=tau3g[iGrid];
    KGLTemplateMetric(Sn,npoint,tau0,tau3,fmin,fmax,fs,vr,gml);
    KGLTemplateMetricEigen(gml,&theta,&lambda1,&lambda2);
    KGLReprodHexagonGrids(1,mmeff,lambda1,lambda2,theta,tau0,tau3,tau0d,tau3d);
    tau0=tau0d[2];
    tau3=tau3d[2];
    KGLTau03ToM12(tau0,tau3,fmin,&m1,&m2);
    if(m1>m1min && m1<m1max && m2>m2min && m2<m2max){
      nGrid2++;
      tau0g[nGrid2]=tau0;
      tau3g[nGrid2]=tau3;
      if(nGrid2 % 500 ==0){
	printf("nGrid=%6d tau0=%7.2f tau3=%5.2f\n",nGrid2,tau0,tau3);
      }
    }
  }
  
  // Then, the largest mass region
  tau0=tau0gmin;
  nGrid2++;
  tau0g[nGrid2]=tau0;
  tau3g[nGrid2]=4.*A3*pow(tau0/4./A0,2/5.);
  KGLTau03ToMcEta(tau0,tau3,fmin,&Mc,&eta);
  Mt=Mc*pow(eta,-3/5.);
  while(Mt<Mtmax){
    KGLTemplateMetric(Sn,npoint,tau0,tau3,fmin,fmax,fs,vr,gml);
    
// dtau0=r03*dtau3 on eta=1/4.
    r03=5/2.*pow(tau3,3/2.)*4*A0/pow(4*A3,5/2.);
    dtau3=(1-mm)/(gml[1][1]*pow(r03,2)+2.*gml[1][2]*r03+gml[2][2]);
    dtau3=-sqrt(fabs(dtau3)); // minus sign: to decrease tau0
    dtau0=r03*dtau3;
    
    //    tau0=tau0+2.*dtau0;
    tau0=tau0+1.5*dtau0;
    tau3=4.*A3*pow(tau0/4./A0,2/5.);
    KGLTau03ToMcEta(tau0,tau3,fmin,&Mc,&eta);
    Mt=Mc*pow(eta,-3/5.);
    
    nGrid2++;
    tau0g[nGrid2]=tau0;
    tau3g[nGrid2]=tau3;
  }
  printf("nGrid=%6d tau0=%7.2f tau3=%5.2f\n",nGrid2,tau0,tau3);
  
  KGLTimeStamp(clock0);
  
  printf("Then, write the grid info...\n");
  char *fname_out1=argv[1];
  if((fout = fopen(fname_out1,"w")) == NULL) {
    perror(fname_out1);
    abort();
  }
  for(iGrid=1;iGrid<=nGrid2;iGrid++){
    fprintf(fout,"%11.6f %9.6f\n",tau0g[iGrid],tau3g[iGrid]);
  }
  fclose(fout);

  char *fname_out2=argv[2];
  if((fout = fopen(fname_out2,"w")) == NULL) {
    perror(fname_out2);
    abort();
  }
  for(iGrid=1;iGrid<=nGrid;iGrid++){
    for(iVertex=1;iVertex<=nVertexMax;iVertex++){
      fprintf(fout,"%f %f\n",tau0v1[iGrid][iVertex],tau3v1[iGrid][iVertex]);
    }
  }
  fclose(fout);
  
  KGLTimeStamp(clock0);
  printf("Then, write the ellipse info...\n");
  
  // Draw isomatch contour (ellipse)
  char *fname_out3=argv[3];
  if((fout = fopen(fname_out3,"w")) == NULL) {
    perror(fname_out3);
    abort();
  }
  int ibin,nbin=100;
  double xwid,xbin;
  for(iGrid=1;iGrid<=nGrid2;iGrid++){
    if(tau0g[iGrid]>0.01*tau0max){
      mmeff=mm;
    }else{
      mmeff=mm*1.01;
    }
    if(tau0g[iGrid]<10. || tau0g[iGrid]>1350.){
      //    if(tau0g[ID]>900. && tau0g[ID]<1100.){
      //    if(tau0g[iGrid]>195.){
      KGLTemplateMetric(Sn,npoint,tau0g[iGrid],tau3g[iGrid],fmin,fmax,fs,vr,gml);
      xwid=gml[2][2]*(1.-mmeff)/(gml[1][1]*gml[2][2]-gml[1][2]*gml[2][1]);
      xwid=2.*sqrt(xwid);
      xbin=xwid/nbin;
      for(ibin=0;ibin<=nbin;ibin++){
        dtau0=-xwid/2.+ibin*xbin;
        dtau3p=pow(gml[1][2]*dtau0,2)-gml[2][2]*(gml[1][1]*pow(dtau0,2)-(1.-mmeff));
        if(dtau3p<0) dtau3p=0;
        
        // upper half
        dtau3=(-gml[1][2]*dtau0+sqrt(dtau3p))/gml[2][2];
        fprintf(fout,"%f %f\n",tau0g[iGrid]+dtau0,tau3g[iGrid]+dtau3);
        
        // lower half
        dtau3=(-gml[1][2]*dtau0-sqrt(dtau3p))/gml[2][2];
        fprintf(fout,"%f %f\n",tau0g[iGrid]+dtau0,tau3g[iGrid]+dtau3);
      }
    }
  }
  fclose(fout);
  
  return 0;
}
