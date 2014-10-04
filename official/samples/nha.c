/* Iterative Least Square fit of given time series */

#include <kagali/KGLStdlib.h>
#include <kagali/nrutil.h>
#include <kagali/KGLTemplateWaveform.h>
#include <kagali/KGLParameters.h>
#include <kagali/KGLIterativeLeastSquare3DNewton.h>
#include <kagali/KGLLeastSquareMethod.h>
#include <kagali/KGLLeastSquareFunc.h>
#include <time.h>
#include <kagali/KGLTimeStamp.h>


int main(int argc,char *argv[])
{
  int i,k;
  int npoint,ishift,nframe,nstart,nend,nshift,kshift,nitr;
  int nfile;
  int nsig; 
  double T,tmin,tmax;
  double Fcostthr;
  double mu0,nu0;
  double fs,df_frame;
  double a2thr;
  
  FILE *fp=NULL;
  clock_t clock0;
  
  if(argc!=16){
    fprintf(stderr,"argc=%d\n",argc);
    fprintf(stderr,"Usage: ./main.out waveform.dat result.dat "
	    "fs tmin tmax nframe nstart nend nshift Fthr a2thr nsig nitr mu0 nu0\n");
    exit(EXIT_FAILURE);
  }
  
  nfile = 2;
  fs=atof(argv[nfile+1]);
  tmin=atof(argv[nfile+2]);
  tmax=atof(argv[nfile+3]);
  T=tmax-tmin;
  nframe=atoi(argv[nfile+4]);
  nstart=atoi(argv[nfile+5]);
  nend=atoi(argv[nfile+6]);
  nshift=atoi(argv[nfile+7]);
  Fcostthr=atof(argv[nfile+8]);
  Fcostthr=pow(10,Fcostthr);
  a2thr=atof(argv[nfile+9]);
  nsig=atoi(argv[nfile+10]);
  nitr=atoi(argv[nfile+11]);
  mu0=atof(argv[nfile+12]);
  nu0=atof(argv[nfile+13]);  
  npoint=(int)(fs*T);
  df_frame=fs/nframe;
  
  printf("You set...\n");
  printf("fs      : %f\n",fs);
  printf("df      : %f\n",df_frame);
  printf("T       : %f\n",T);
  printf("npoint  : %d\n",npoint);
  printf("nframe  : %d\n",nframe);
  printf("nstart  : %d\n",nstart);
  printf("nend    : %d\n",nend);
  printf("nshift  : %d\n",nshift);
  printf("Fcostthr: %e\n",Fcostthr);
  printf("a2thr   : %e\n",a2thr);
  printf("nsig    : %d\n",nsig);
  printf("nitr    : %d\n",nitr);
  printf("mu0     : %f\n",mu0);
  printf("nu0     : %f\n",nu0);
  
  clock0 = clock();
  KGLTimeStamp(clock0);
  
  double *datat;
  datat=dvector(0,npoint-1);
  char *fname1=argv[1];
  fp=fopen(fname1,"rb");
  for(i=0;i<npoint;i++){
    fread(&datat[i], sizeof(double), 1, fp);
    //    printf("%e %e\n",time[i],datat[i]);
  }
  fclose(fp);
  
  double *frame;
  double *Afit,*pfit,*ffit;
  frame=dvector(0,nframe-1);
  Afit=dvector(0,nsig-1);
  pfit=dvector(0,nsig-1);
  ffit=dvector(0,nsig-1);
  
  char *fname2=argv[2];
  fp=fopen(fname2,"w");
  ishift=0;
  for(kshift=nstart;kshift<nend;kshift+=nshift){
    if(kshift+nframe>npoint){
      printf("kshift + nframe > npoint\n");
      break;
    }
    if(kshift%100 == 0) printf("kshift=%d\n",kshift);
    for(k=0;k<nframe;k++){
      frame[k]=datat[k+kshift];
    }
    KGLIterativeLeastSquare3DNewton(frame,fs,nframe,Fcostthr,a2thr,nsig,nitr,mu0,nu0,Afit,pfit,ffit);
    for(i=0;i<nsig;i++){
      if(Afit[i]!=0){
	printf("%e %e %e %e %e\n",kshift/fs,(kshift+nframe)/fs,Afit[i],pfit[i],ffit[i]);
	fprintf(fp,"%e %e %e %e\n",tmin+(kshift+nframe/2.)/fs,Afit[i],pfit[i],ffit[i]);
      }
    }
  }
  fclose(fp);
  
  return 0;
}
