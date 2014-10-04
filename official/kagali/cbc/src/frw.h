int fwrite_svd(char fname[],int ngrid,int nbasis,double srate,int npoint,
	       int *nth,double *time,double *normh,
	       double *s,double **v,double **u);

int fwrite_htemplate(char fname[],int ngrid,int npoint,
		     double *time,double *normh,double **htemplate);

int fwrite_para(char fname[],int ngrid,double *tau0g,double *tau3g,
		double *Mtotg,double *etag,double *m1g,double *m2g);

int fread_slice(char fname[],int *ngrid,int *nbasis,double *srate,int *npoint);

int fread_svd(char fname[],int *ngrid,int *nbasis,double *srate,int *npoint,
	      int *nth,double *time,double *normh,
	      double *s,double **v,double **u);

int fread_htemplate(char fname[],int ngrid,int npoint,
		    double *time,double *normh,double **htemplate);

