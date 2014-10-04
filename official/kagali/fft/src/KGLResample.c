/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *  This program changes sampling rate by rational factor.
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
#include <time.h>
#include <fcntl.h> 
#include <kagali/nrutil.h>
#include <kagali/KGLBandPassFilterUeno.h>
#include <kagali/KGLParameters.h>


int KGLResample(double *time,double *ht,double *timed,double *htd,
		double fs,double fs3,int npoint,int iwindow,int nkernel)
{
  int j,k;
  int npointi,npointd;
  double fs2;
  int irrate1,irrate2;
  double rrate1,rrate2;
  int gcd,lcm;
  int ifs,ifs3;
  int ifstmp,ifs3tmp;
  int itmp;
  ifs=(int)(fs);
  ifs3=(int)(fs3);
  ifstmp=ifs;
  ifs3tmp=ifs3;
  while(ifs3tmp != 0){
    itmp=ifs3tmp;
    ifs3tmp=ifstmp%ifs3tmp;
    ifstmp=itmp;
  }
  gcd=ifstmp;
  lcm=(ifs*ifs3)/gcd;
  fs2=(double)lcm;
  rrate1=fs2/fs;
  rrate2=fs2/fs3;
  irrate1=(int)(rrate1);
  irrate2=(int)(rrate2);
  npointi=npoint*irrate1;
  npointd=(int)(fs3/fs*(double)npoint);  
  
  double tmp1;
  double *timei,*htz;
  timei=dvector(0,npointi-1);
  htz=dvector(0,npointi-1);  
  k=0;
  for(tmp1=time[0];tmp1<=time[npoint-1];tmp1+=1/fs2){
    timei[k]=tmp1;
    if(k%irrate1==0){
      htz[k]=ht[k/irrate1];
    }else{
      htz[k]=0;
    }
    k++;
  }
  
  double *lpfh;
  lpfh=dvector(0,nkernel-1);  
  if(ifs==ifs3){
    for(k=0;k<npointd;k++){
      timed[k]=timei[k];
      htd[k]=htz[k];
    }
  }else if(ifs!=ifs3){
    if(ifs<ifs3){
      lpf_norm(iwindow,nkernel,fs/fs2/2.,lpfh);
    }else if(ifs>ifs3){
      lpf_norm(iwindow,nkernel,fs3/fs2/2.,lpfh);
    }
    for(k=0;k<npointd;k++){
      timed[k]=timei[k*irrate2];
      for(j=0;j<nkernel;j++){
	if(k*irrate2+(nkernel-1)/2-j>=0 && k*irrate2+(nkernel-1)/2-j<npointi){
	  htd[k]+=rrate1*htz[k*irrate2+(nkernel-1)/2-j]*lpfh[j];
	}
      }
    }
  }
  
  return 0;
}
