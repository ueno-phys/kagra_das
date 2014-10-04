#ifndef TEMPLATE_OPERATION_H
#define TEMPLATE_OPERATION_H

#include <complex.h>

KGL_BEGIN_DECLS

void KGLInnerProdFreq( 
    double complex *a1,
    double complex *a2,
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double complex *out /** only real part should be used*/
    );

void KGLInnerProdTime( 
   double *a,
   double *b,
   int n,
   double *sum
   );

void KGLWhitening( 
    double complex *a,
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double complex *aw
    );


KGL_END_DECLS

#endif /* TEMPLATE_OPERATION_H */
