/**
 * \author K. Oohara
 * \file
 *
 */

#include "KGLVersion.h"
#include <kagali/KGLStdlib.h>
#include <kagali/BandPassFilter.h>

/**
 * simple band-pass filter
 *
 * parameters
 *  status: int; return status information
 *  data  : array of double; gives an original signal
 *                           and returns the filtered signal
 *  dataLen: int; size of array data
 *  dt     : double; sampling interval = 1/(Nyquist freq.)
 *  freqLow : minimum freq.
 *  freqHigh: maximum freq.
 *
 *  This will be low-pass filter, if freqLow <= 0.0.
 *  This will be high-pass filter, if freqHigh <= 0.0.
 */

void KGLBandPassFilter ( //begin{proto}
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double       *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,    /**< (in) > 0                         */
    const double dt,      /**< (in) > 0.0                       */
    const double freqLow, /**< (in) minimum freq. of the filter */
    const double freqHigh /**< (in) maximum freq. of the filter */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,data != NULL,"data is NULL");
    KGLAssert(status,size > 0,"size must be larger than 0");
    KGLAssert(status,dt > 0,"dt must be larger than 0");
    KGLAssert(status,freqLow >= 0,"freqLow is less than 0");
    KGLAssert(status,(freqHigh == 0 || freqHigh > freqLow),
              "freqHigh is less than freqLow");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

    int k, iMin, iMax;
    bool fftonce = false;
    bool fftnew  = false;

    KGLRealFFTPlan *plan1  = NULL;
    double         *rdataP = NULL;
    fftw_complex   *cdataP = NULL;
    fftw_plan      r2c     = NULL;
    fftw_plan      c2r     = NULL;

    if(plan == NULL) {
        /* calculate FFT only once */
        fftonce = true;
        fftnew = true;
    }
    else if(*plan == NULL) {
        /* the first call */
        fftnew = true;
    }

    iMin = (int)(freqLow*dt*size);
    iMax = (int)(freqHigh*dt*size);

    if(iMin > size/2) {
        /* data[*] = 0.0 if freqLow > (Nyquist freq.) */
        for(k = 0; k < size; k++) data[k] = 0.0;
        return;
    }
    if((iMin == 0) && ((iMax == 0) || (iMax >= size/2))) {
        /* nothing to do, if low-pass filter && freqHigh >= (Nyquist freq.) */
        return;
    }

    if(fftnew) {
        plan1 = KGLCreateRealFFTPlan(status,KGL_FFT_BOTH,size);
        if(plan1 == NULL) {
            KGLAddMessage(status,"");
            return;
        }
#ifdef KGL_PTHREAD
        rdataP = (double *)fftw_malloc(sizeof(double)*size);
        cdataP = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size/2+1));
        if((rdataP == NULL) || (cdataP == NULL)) {
            KGLAddError(status,"fail to alloc rdata and/or cdata\n");
            return;
        }
#else
        rdataP  = plan1->rdata;
        cdataP  = plan1->cdata;
#endif
        r2c     = plan1->r2c;
        c2r     = plan1->c2r;
        if(!fftonce) *plan = plan1;
    }
    else {
        plan1   = *plan;
        KGLAssert(status,plan1->size == size,"size != plan->size");
        if(KGLCheckError(status)) return;

        r2c     = plan1->r2c;
        c2r     = plan1->c2r;
        if((r2c == NULL) || (c2r == NULL)) {
            KGLAddRealFFTPlan(status,plan1);
            if(KGLCheckError(status)) return;
            r2c = plan1->r2c;
            c2r = plan1->c2r;
        }

#ifdef KGL_PTHREAD
        rdataP = (double *)fftw_malloc(sizeof(double)*size);
        cdataP = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size/2+1));
        if((rdataP == NULL) || (cdataP == NULL)) {
            KGLAddError(status,"fail to alloc rdata and/or cdata\n");
            return;
        }
#else
        rdataP  = plan1->rdata;
        cdataP  = plan1->cdata;
        KGLAssert(status,(rdataP != NULL) && (cdataP != NULL),
            "plan->rdata == NULL and/or plan->cdata == NULL");
        if(KGLCheckError(status)) return;
#endif
    }

    for(k = 0; k < size; k++) rdataP[k] = data[k];
#ifdef KGL_PTHREAD
    fftw_execute_r2c(r2c,rdataP,cdataP);
#else
    fftw_execute(r2c);
#endif

    if(iMin > 0) {
        for(k = 0; k < iMin; k++) cdataP[k] = 0.0;
    }
    if(iMax > 0 && iMax < size/2) {
        for(k = iMax+1; k <= size/2; k++) cdataP[k] = 0.0;
    }
#ifdef KGL_PTHREAD
    fftw_execute_c2r(c2r,cdataP,rdataP);
#else
    fftw_execute(c2r);
#endif

    for(k = 0; k < size; k++) data[k] = rdataP[k]/(double)size;

    if(fftonce) KGLDestroyRealFFTPlan(plan1);

    return;
}

void KGLLowPassFilter ( //begin{proto}
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,      /**< (in) > 0                         */
    const double dt,        /**< (in) > 0.0                       */
    const double freqHigh   /**< (in) maximum freq. of the filter */
    ) //end{proto}
{
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,freqHigh > 0.0,"freqHigh must be larger than 0");
    if(KGLCheckError(status)) return;

    if(freqHigh >= 0.5/dt) return;

    KGLBandPassFilter(status,plan,data,size,dt,0.0,freqHigh);

    return;
}

void KGLHighPassFilter ( //begin{proto}
    KGLStatus      *status, /**< status                           */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *data,   /**< (in/out) signal: array [dataLen] */
    const size_t size,      /**< (in) > 0                         */
    const double dt,        /**< (in) > 0.0                       */
    const double freqLow    /**< (in) minimum freq. of the filter */
    ) //end{proto}
{
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,freqLow > 0.0,"freqLow must be larger than 0");
    if(KGLCheckError(status)) return;

    if(freqLow > 0.5/dt) {
        for(int k = 0; k < size; k++) data[k] = 0.0;
        return;
    }

    KGLBandPassFilter(status,plan,data,size,dt,freqLow,0.0);

    return;
}
