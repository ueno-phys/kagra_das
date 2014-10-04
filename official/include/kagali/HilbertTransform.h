#ifndef HILBERTTRANSFORM_H
#define HILBERTTRANSFORM_H

#include <kagali/RealFFT.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_spline.h>

KGL_BEGIN_DECLS

void KGLHilbertTransform( 
    KGLStatus      *status, /**< status                 */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *idata,  /**< (out) returns the Hilbert transform */
    const double   *rdata,  /**< (in)  gives a signal */
    const size_t   size     /**< length of data     */
    );

void KGLHilbertTransformOnce( 
    KGLStatus    *status, /**< status                 */
    double       *idata,  /**< (out) returns the Hilbert transform */
    const double *rdata,  /**< (in)  gives a signal */
    const size_t size     /**< length of data     */
    );

void KGLHSAFreqAmp( 
    KGLStatus    *status, /**< status                 */
    double       *phi,    /**< (out): returns the phase */
    double       *freq,   /**< (out): returns the instantaneous frequency */
    double       *amp,    /**< (out): returns the instantaneous amplitude */
    const double *t,      /**< (in): gives the time */
    const double *rdata,  /**< (in): gives the real part of data */
    const double *idata,  /**< (in): gives the imaginary part of data */
    const int size        /**< length of data         */
    );

void KGLHSAFrequency( 
    KGLStatus    *status, /**< status                 */
    double       *freq, /**< (out): array [size] */
    const double *t,    /**< (in): array [size]  */
    const double *phi,  /**< (in): array [size]  */
    const size_t size,  /**< length of data         */
    const int nderiv    /**< the number of points for calculating derivatives */
    );


KGL_END_DECLS

#endif /* HILBERTTRANSFORM_H */
