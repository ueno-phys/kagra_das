#ifndef ILS_H
#define ILS_H

KGL_BEGIN_DECLS

void KGLIterativeLeastSquare3DNewton( 
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
   );


KGL_END_DECLS

#endif /* ILS_H */
