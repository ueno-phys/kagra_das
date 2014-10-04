#include <kagali/KGLStdlib.h>
#include <kagali/nrutil.h>
#include <kagali/KGLRWBinary.h>

int KGLFwriteSVD( //begin{proto}
   char fname[],
   int ngrid,
   int nbasis,
   double srate,
   int npoint,
   int *nth,
   double *time,
   double *normh,
   double *s,
   double **v,
   double **u
   ) //end{proto}
{
  int i,j;
  FILE *fp;
  fp=fopen(fname,"wb");
  fwrite(&ngrid, sizeof(int), 1, fp);
  fwrite(&nbasis, sizeof(int), 1, fp);
  fwrite(&srate, sizeof(double), 1, fp);
  fwrite(&npoint, sizeof(int), 1, fp);
  for(i=1;i<=6;i++){
    fwrite(&nth[i], sizeof(int), 1, fp);
  }
  for(j=0;j<npoint;j++){
    fwrite(&time[j], sizeof(double), 1, fp);
  }
  for(i=1;i<=2*ngrid;i++){
    fwrite(&normh[i], sizeof(double), 1, fp);
    fwrite(&s[i], sizeof(double), 1, fp);
    for(j=1;j<=2*ngrid;j++){
      fwrite(&v[i][j], sizeof(double), 1, fp);
    }
    if(i<=nbasis) {
      for(j=0;j<npoint;j++){
	fwrite(&u[i][j], sizeof(double), 1, fp);
      }
    }
  }
  fclose(fp);
  
  return 0;
}


int KGLFwriteHtemplate( //begin{proto}
   char fname[],
   int ngrid,
   int npoint,
   double *time,
   double *normh,
   double **htemplate
   ) //end{proto}
{
  int i,j;
  FILE *fp;
  fp=fopen(fname,"wb");
  for(j=0;j<npoint;j++){
    fwrite(&time[j], sizeof(double), 1, fp);
  }
  for(i=1;i<=2*ngrid;i++){
    fwrite(&normh[i], sizeof(double), 1, fp);
    for(j=0;j<npoint;j++){
      fwrite(&htemplate[i][j], sizeof(double), 1, fp);
    }
  }
  fclose(fp);
  
  return 0;
}


int KGLFwritePara( //begin{proto}
   char fname[],
   int ngrid,
   double *tau0g,
   double *tau3g,
   double *Mtotg,
   double *etag,
   double *m1g,
   double *m2g
   ) //end{proto}
{
  int i;
  FILE *fp;
  fp=fopen(fname,"wb");
  for(i=1;i<=ngrid;i++){
    fwrite(&tau0g[i], sizeof(double), 1, fp);
    fwrite(&tau3g[i], sizeof(double), 1, fp);
    fwrite(&Mtotg[i], sizeof(double), 1, fp);
    fwrite(&etag[i], sizeof(double), 1, fp);
    fwrite(&m1g[i], sizeof(double), 1, fp);
    fwrite(&m2g[i], sizeof(double), 1, fp);
  }
  fclose(fp);
  
  return 0;
}


int KGLFreadSlice( //begin{proto}
   char fname[],
   int *ngrid_tmp,
   int *nbasis_tmp,
   double *srate_tmp,
   int *npoint_tmp
   ) //end{proto}
{
  int ngrid,nbasis,npoint;
  double srate;
  FILE *fp;
  if((fp = fopen(fname,"rb"))==NULL) {
    printf("open error.\n");
    exit(1);
  }
  fread(&ngrid, sizeof(int), 1, fp);
  fread(&nbasis, sizeof(int), 1, fp);
  fread(&srate, sizeof(double), 1, fp);
  fread(&npoint, sizeof(int), 1, fp);
  *ngrid_tmp=ngrid;
  *nbasis_tmp=nbasis;
  *srate_tmp=srate;
  *npoint_tmp=npoint;
  
  fclose(fp);
  
  return 0;
}


int KGLFreadSVD( //begin{proto}
   char fname[],
   int *ngrid_tmp,
   int *nbasis_tmp,
   double *srate_tmp,
   int *npoint_tmp,
   int *nth,
   double *time,
   double *normh,
   double *s,
   double **v,
   double **u
   )  //end{proto}
{
  int i,j;
  int ngrid,nbasis,npoint;
  double srate;
  FILE *fp;
  if((fp = fopen(fname,"rb"))==NULL) {
    printf("open error.\n");
    exit(1); 
 }
  fread(&ngrid, sizeof(int), 1, fp);
  fread(&nbasis, sizeof(int), 1, fp);
  fread(&srate, sizeof(double), 1, fp);
  fread(&npoint, sizeof(int), 1, fp);
  *ngrid_tmp=ngrid;
  *nbasis_tmp=nbasis;
  *srate_tmp=srate;
  *npoint_tmp=npoint;
  for(i=1;i<=6;i++){
    fread(&nth[i], sizeof(int), 1, fp);
  }
  for(j=0;j<npoint;j++){
    fread(&time[j], sizeof(double), 1, fp);
  }
  for(i=1;i<=2*ngrid;i++){
    fread(&normh[i], sizeof(double), 1, fp);
    fread(&s[i], sizeof(double), 1, fp);
    for(j=1;j<=2*ngrid;j++){
      fread(&v[i][j], sizeof(double), 1, fp);
    }
    if(i<=nbasis) {
      for(j=0;j<npoint;j++){
	fread(&u[i][j], sizeof(double), 1, fp);
      }
    }
  }
  fclose(fp);
  
  return 0;
}


int KGLFreadHtemplate( //begin{proto}
   char fname[],
   int ngrid,
   int npoint,
   double *time,
   double *normh,
   double **htemplate
   ) //end{proto}
{
  int i,j;
  FILE *fp;
  if((fp = fopen(fname,"rb"))==NULL) {
    printf("open error.\n");
    exit(1); 
  }
  time=dvector(0,npoint-1);
  for(j=0;j<npoint;j++){
    fread(&time[j], sizeof(double), 1, fp);
  }
  for(i=1;i<=2*ngrid;i++){
    fread(&normh[i], sizeof(double), 1, fp);
    for(j=0;j<npoint;j++){
      fread(&htemplate[i][j], sizeof(double), 1, fp);
    }
  }
  fclose(fp);
  
  return 0;
}
