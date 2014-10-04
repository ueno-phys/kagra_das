#ifndef BANDPASSFILTER_H
#define BANDPASSFILTER_H

#include <kagali/RealFFT.h>
#include <complex.h>
#include <fftw3.h>

KGL_BEGIN_DECLS

void KGLBandPassFilter ( 
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double       *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,    /**< (in) > 0                         */
    const double dt,      /**< (in) > 0.0                       */
    const double freqLow, /**< (in) minimum freq. of the filter */
    const double freqHigh /**< (in) maximum freq. of the filter */
    );

void KGLLowPassFilter ( 
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,      /**< (in) > 0                         */
    const double dt,        /**< (in) > 0.0                       */
    const double freqHigh   /**< (in) maximum freq. of the filter */
    );

void KGLHighPassFilter ( 
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,      /**< (in) > 0                         */
    const double dt,        /**< (in) > 0.0                       */
    const double freqLow    /**< (in) minimum freq. of the filter */
    );


KGL_END_DECLS

#endif /* BANDPASSFILTER_H */
