#ifndef RW_BINAY_H
#define RW_BINAY_H

KGL_BEGIN_DECLS

int KGLFwriteSVD( 
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
   );

int KGLFwriteHtemplate( 
   char fname[],
   int ngrid,
   int npoint,
   double *time,
   double *normh,
   double **htemplate
   );

int KGLFwritePara( 
   char fname[],
   int ngrid,
   double *tau0g,
   double *tau3g,
   double *Mtotg,
   double *etag,
   double *m1g,
   double *m2g
   );

int KGLFreadSlice( 
   char fname[],
   int *ngrid_tmp,
   int *nbasis_tmp,
   double *srate_tmp,
   int *npoint_tmp
   );

int KGLFreadSVD( 
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
   );

int KGLFreadHtemplate( 
   char fname[],
   int ngrid,
   int npoint,
   double *time,
   double *normh,
   double **htemplate
   );


KGL_END_DECLS

#endif /* RW_BINAY_H */
