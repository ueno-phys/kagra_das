/**
 * \author K. Oohara
 * \file
 *
 * \brief Functions to peform the Hilbert transform
 *
 * 
 */

#include "KGLVersion.h"
#include <kagali/KGLStdlib.h>
#include <kagali/HilbertTransform.h>
#include <string.h>

/**
 * to calculate Hilbert transform
 *
 * HTPlan **plan: 
 *    - *plan must be NULL on the first call. The created plan will be
 *      stored in *plan.
 *    - If *plan != NULL, this function will use fftw_plan and data arrays
 *      which are kept in *plan.
 *    - If plan == NULL, this function will not keep the plan.
 *    - The plan should be freed by the user program 
 *      using XKGLDestroyRealHTPlan() unless plan == NULL.
 *
 */
void KGLHilbertTransform( //begin{proto}
    KGLStatus      *status, /**< status                 */
    KGLRealFFTPlan **plan,  /**< (in/out) FFTPlan will be stored in *plan */
    double         *idata,  /**< (out) returns the Hilbert transform */
    const double   *rdata,  /**< (in)  gives a signal */
    const size_t   size     /**< length of data     */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,idata != NULL,"idata is NULL");
    KGLAssert(status,rdata != NULL,"rdata is NULL");
    KGLAssert(status,size > 0,"size must be larger than 0");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

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
        if(plan1->size != size) {
            KGLAddError(status,"plan->size != size\n");
            return;
        }

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
        if((rdataP == NULL) || (cdataP == NULL)) {
            KGLAddError(status,"(plan->[rc]data == NULL)\n");
            return;
        }
#endif
    }


//    for(int k = 0; k < size; k++) rdataP[k] = rdata[k];
    memcpy(rdataP,rdata,sizeof(double)*size);
    fftw_execute(r2c);

    cdataP[0] = 0;
    cdataP[size/2] = 0;

    for(int k = 1; k < size/2; k++) cdataP[k] *= -I/(double)size;
    fftw_execute(c2r);
//    for(int k = 0; k < size; k++) idata[k] = rdataP[k];
    memcpy(idata,rdataP,sizeof(double)*size);

    if(fftonce) KGLDestroyRealFFTPlan(plan1);

    return;
}

void KGLHilbertTransformOnce( //begin{proto}
    KGLStatus    *status, /**< status                 */
    double       *idata,  /**< (out) returns the Hilbert transform */
    const double *rdata,  /**< (in)  gives a signal */
    const size_t size     /**< length of data     */
    ) //end{proto}
{
    KGLHilbertTransform(status,NULL,idata,rdata,size);
    return;
}

/**
 *  Calculate the phase and the instantaneous frequency and amplitude.
 *  phi[i] = atan2(datai[i],datar[i]) + 2*pi*n
 *    the integer n is chosen so that phi will be continuous.
 *  freq[i] = d(phi[i])/dt
 *          = (datar[i]*d(datai[i])/dt - datai[i]*d(datar[i])/dt)
 *            /(datar[i]*datar[i] + datai[i]*datai[i])
 */

void KGLHSAFreqAmp( //begin{proto}
    KGLStatus    *status, /**< status                 */
    double       *phi,    /**< (out): returns the phase */
    double       *freq,   /**< (out): returns the instantaneous frequency */
    double       *amp,    /**< (out): returns the instantaneous amplitude */
    const double *t,      /**< (in): gives the time */
    const double *rdata,  /**< (in): gives the real part of data */
    const double *idata,  /**< (in): gives the imaginary part of data */
    const int size        /**< length of data         */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,phi != NULL,"phi is NULL");
    KGLAssert(status,freq != NULL,"freq is NULL");
    KGLAssert(status,t != NULL,"t is NULL");
    KGLAssert(status,rdata != NULL,"rdata is NULL");
    KGLAssert(status,idata != NULL,"idata is NULL");
    KGLAssert(status,size > 1,"size must be larger than 1");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

    int it;
    double phi0, phip, phii;

    it = 0;
    phi0 = 0.0;
    phip = atan2(idata[it],rdata[it]);
    phi[it] = phip;
    for(it = 1; it < size; it++) {
        phii = atan2(idata[it],rdata[it]);
        if(phii+phi0-phip < -KGL_PI) {
            phi0 += 2*KGL_PI;
        }
        else if(phii+phi0-phip > KGL_PI) {
            phi0 -= 2*KGL_PI;
        }
        phip = phii + phi0;
        phi[it] = phip;
    }

    gsl_interp *interp = gsl_interp_alloc(gsl_interp_cspline, size);
    if(interp == NULL) {
        KGLAddError(status,"fail to alloc gsl_interp_alloc\n");
        return;
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    if(acc == NULL) {
        KGLAddError(status,"fail to alloc gsl_interp_accel");
        return;
    }
    gsl_interp_init(interp,t,phi,size);
    for(it = 0; it < size; it++) {
        freq[it] = gsl_interp_eval_deriv(interp,t,phi,t[it],acc)/(2*KGL_PI);
    }
    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);

    for(it = 0; it < size; it++) {
        amp[it] = sqrt(SQ(rdata[it]) + SQ(idata[it]));
    }

    return;
}

void KGLHSAFrequency( //begin{proto}
    KGLStatus    *status, /**< status                 */
    double       *freq, /**< (out): array [size] */
    const double *t,    /**< (in): array [size]  */
    const double *phi,  /**< (in): array [size]  */
    const size_t size,  /**< length of data         */
    const int nderiv    /**< the number of points for calculating derivatives */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,t != NULL,"t is NULL");
    KGLAssert(status,phi != NULL,"phi is NULL");
    KGLAssert(status,size > 1,"size must be larger than 1");
    KGLAssert(status,nderiv > 1,"nderiv must be larger than 1");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

    int idr, it, it0, it1;
    double dt;
    double *coeff;

    coeff = (double *)calloc(nderiv,sizeof(double));
    if(coeff == NULL) {
        KGLAddError(status,"fail to alloc coeff\n");
        return;
    }

    dt = t[1] - t[0];
    for(idr = 0; idr < nderiv; idr++) {
        coeff[idr] = 6.0*(2.0*(double)(idr+1)-(double)(nderiv+1))
            /(dt*(double)nderiv*(double)(nderiv-1)*(double)(nderiv+1));
    }

    it0 = nderiv-2-(nderiv-1)/2;
    it1 = size-(nderiv-1)+it0;
    for(it = it0; it < it1; it++) {
        double omega = 0.0;
        for(idr = 0; idr < nderiv; idr++) {
            omega += coeff[idr]*phi[it+idr-it0];
        }
        freq[it] = omega/(2*KGL_PI);
    }
    for(it = 0; it < it0; it++) freq[it] = freq[it0];
    for(it = it1; it < size; it++) freq[it] = freq[it1-1];

    return;
}
