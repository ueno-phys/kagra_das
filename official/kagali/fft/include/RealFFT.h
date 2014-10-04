#ifndef REALFFT_H
#define REALFFT_H

#include <complex.h>
#include <fftw3.h>

KGL_BEGIN_DECLS

#define KGL_FFT_FORWARD (1 << 0)
#define KGL_FFT_REVERSE (1 << 1)
#define KGL_FFT_BOTH (KGL_FFT_FORWARD|KGL_FFT_REVERSE)

struct KGLRealFFTPlan_tag {
    size_t       size;     /**< length of the real data vector for this plan */
    double       *rdata;   /**< array of the real data */
    fftw_complex *cdata;   /**< array of the complex data */
    fftw_plan    r2c;      /**< FFTW plan for forward transform */
    fftw_plan    c2r;      /**< FFTW plan for reverse transform */
};

typedef struct KGLRealFFTPlan_tag KGLRealFFTPlan;

KGLRealFFTPlan *KGLCreateRealFFTPlan( 
    KGLStatus *status,  /**< status */
    const int direct,   /**< direction flag */
    const size_t size   /**< length of the real data vector */
    );

void KGLAddRealFFTPlan( 
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan *plan     /**< (in/out) KGLRealFFTPlan to be changed */
    );

void KGLDestroyRealFFTPlan(KGLRealFFTPlan *plan);

void KGLForwardRealFFT( 
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan **plan,   /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double complex *cdata,   /**< (out) returns the Fourier transform */
    const double   *rdata,   /**< (in) gives a real data */
    const size_t   size      /**< (in) gives the size of the real data */
    );

void KGLForwardRealFFTOnce( 
    KGLStatus      *status, /**< status */
    double complex *cdata,  /**< (out) returns the Fourier transform */
    const double   *rdata,  /**< (in) gives a real data */
    const size_t   size     /**< (in) gives the size of the real data */
    );

void KGLReverseRealFFT( 
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan **plan,   /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double         *rdata,   /**< (out) returns the reverse Fourier transform */
    const double complex *cdata, /**< (in) gives a complex data */
    const size_t   size      /**< (in) gives the size of the real data */
    );

void KGLReverseRealFFTOnce( 
    KGLStatus      *status,  /**< status */
    double         *rdata,   /**< (out) returns the reverse Fourier transform */
    const double complex *cdata, /**< (in) gives a complex data */
    const size_t   size      /**< (in) gives the size of the real data */
    );


KGL_END_DECLS

#endif /* REALFFT_H */
