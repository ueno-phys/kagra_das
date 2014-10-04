#include "KGLVersion.h"
#include <kagali/KGLStdlib.h>
#include "configRealFFT.h"
#include <kagali/RealFFT.h>
#include <string.h>

static fftw_plan xkgl_fftw_plan_dft_r2c_1d(
    int n0, double *in, fftw_complex *out
    );

static fftw_plan xkgl_fftw_plan_dft_r2c_1d(
    int n0, double *in, fftw_complex *out
    )
{
    uint32_t flag;
    fftw_plan plan;

    KGL_FFTW_LOCK;

    flag = KGL_FFTWPLAN|FFTW_WISDOM_ONLY;
    plan = fftw_plan_dft_r2c_1d(n0,in,out,flag);
    if(plan != NULL) goto ret;

#if 0
    if(fftw_import_system_wisdom() != 0) {
        plan = fftw_plan_dft_r2c_1d(n0,in,out,flag);
        if(plan != NULL) goto ret;
    }
#endif

    flag = KGL_FFTWPLAN_FS;
    plan = fftw_plan_dft_r2c_1d(n0,in,out,flag);

ret:
    KGL_FFTW_UNLOCK;

    return plan;

}

static fftw_plan xkgl_fftw_plan_dft_c2r_1d(
    int n0, fftw_complex *in, double *out
    );

static fftw_plan xkgl_fftw_plan_dft_c2r_1d(
    int n0, fftw_complex *in, double *out
    )
{
    uint32_t flag;
    fftw_plan plan;

    KGL_FFTW_LOCK;

    flag = KGL_FFTWPLAN|FFTW_WISDOM_ONLY;
    plan = fftw_plan_dft_c2r_1d(n0,in,out,flag);
    if(plan != NULL) goto ret;

#if 0
    if(fftw_import_system_wisdom() != 0) {
        plan = fftw_plan_dft_c2r_1d(n0,in,out,flag);
        if(plan != NULL) goto ret;
    }
#endif

    flag = KGL_FFTWPLAN_FS;
    plan = fftw_plan_dft_c2r_1d(n0,in,out,flag);

ret:
    KGL_FFTW_UNLOCK;

    return plan;

}

/**
 * \section sec_RealFFT_KAGALI
 *
 * KAGALI wrapper of the FFTW3 fast Fourier transform
 *
 * Real forward/reverse FFT
 *
 */

/**
 *  to create KGLRealFFTPlan
 *
 *  RETURN VALUE
 *   KGLCreateRealFFTPlan() returns a pointer to KGLRealFFTPlan.
 *   On error, it returns NULL.
 */
KGLRealFFTPlan *KGLCreateRealFFTPlan( //begin{proto}
    KGLStatus *status,  /**< status */
    const int direct,   /**< direction flag */
    const size_t size   /**< length of the real data vector */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,(direct & KGL_FFT_FORWARD) || (direct & KGL_FFT_REVERSE),
              NULL);
    KGLAssert(status,size > 0,"size must be larger than 0");
    if(KGLCheckError(status)) return NULL;
    /* end of check on parameters */

    KGLRealFFTPlan *plan  = NULL;
    double         *rdata = NULL;
    fftw_complex   *cdata = NULL;
    fftw_plan      r2c    = NULL;
    fftw_plan      c2r    = NULL;

    plan = (KGLRealFFTPlan *)calloc(1,sizeof(KGLRealFFTPlan));
    if(plan == NULL) {
        KGLAddError(status,"fail to alloc plan\n");
        return NULL;
    }
    rdata = (double *)fftw_malloc(sizeof(double)*size);
    cdata = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size/2+1));
    if((rdata == NULL) || (cdata == NULL)) {
        KGLAddError(status,"fail to alloc rdata and/or cdata\n");
        KGLDestroyRealFFTPlan(plan);
        return NULL;
    }

    if(direct & KGL_FFT_REVERSE) {
        c2r = xkgl_fftw_plan_dft_c2r_1d(size,cdata,rdata);
        if(c2r == NULL) {
            KGLAddError(status,"fail to fftw_plan_dft_c2r\n");
            KGLDestroyRealFFTPlan(plan);
            return NULL;
        }
    }
    if(direct & KGL_FFT_FORWARD) {
        r2c = xkgl_fftw_plan_dft_r2c_1d(size,rdata,cdata);
        if(r2c == NULL) {
            KGLAddError(status,"fail to fftw_plan_dft_r2c\n");
            KGLDestroyRealFFTPlan(plan);
            return NULL;
        }
    }

#ifdef KGL_PTHREAD
    fftw_free(rdata);
    fftw_free(cdata);
    rdata = NULL;
    cdata = NULL;
#endif

    plan->size   = size;
    plan->rdata  = rdata;
    plan->cdata  = cdata;
    plan->r2c    = r2c;
    plan->c2r    = c2r;

    return plan;
}

void KGLAddRealFFTPlan( //begin{proto}
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan *plan     /**< (in/out) KGLRealFFTPlan to be changed */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,plan != NULL,"plan is NULL");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

    double *rdata;
    fftw_complex *cdata;

    size_t size  = plan->size;
    if(size <= 0) {
        KGLAddError(status,"plan->size <= 0\n");
        return;
    }

#ifdef KGL_PTHREAD
    rdata = (double *)fftw_malloc(sizeof(double)*size);
    cdata = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size/2+1));
    if((rdata == NULL) || (cdata == NULL)) {
        KGLAddError(status,"fail to alloc rdata and/or cdata\n");
        return;
    }
#else
    rdata = plan->rdata;
    cdata = plan->cdata;
    if((rdata == NULL) || (cdata == NULL)) {
        KGLAddError(status,"(plan->[rc]data == NULL)\n");
        return;
    }
#endif

    fftw_plan r2c = plan->r2c;
    fftw_plan c2r = plan->c2r;

    if(r2c == NULL) {
        r2c = xkgl_fftw_plan_dft_r2c_1d(size,rdata,cdata);
        if(r2c == NULL) {
            KGLAddError(status,"fail to fftw_plan_dft_r2c\n");
            return;
        }
        plan->r2c = r2c;
    }

    if(c2r == NULL) {
        c2r = xkgl_fftw_plan_dft_c2r_1d(size,cdata,rdata);
        if(c2r == NULL) {
            KGLAddError(status,"fail to fftw_plan_dft_c2r\n");
            return;
        }
        plan->c2r = c2r;
    }

#ifdef KGL_PTHREAD
    fftw_free(rdata);
    fftw_free(cdata);
#endif

    return;

}


void KGLDestroyRealFFTPlan(KGLRealFFTPlan *plan) //prototype
{
    if(plan == NULL) return;

    KGL_FFTW_LOCK;
    fftw_destroy_plan(plan->r2c);
    fftw_destroy_plan(plan->c2r);
    KGL_FFTW_UNLOCK;

    fftw_free(plan->rdata);
    fftw_free(plan->cdata);

    free(plan);

    return;
}

/**
 * to calculate forward FFT
 *
 * KGLRealFFTPlan **plan: 
 *    - *plan must be NULL on the first call. The created plan will be
 *      stored in *plan.
 *    - If *plan != NULL, this function will use fftw_plan and data arrays
 *      which are kept in *plan.
 *    - If plan == NULL, this function will not keep the plan.
 *    - The plan should be freed by the user program 
 *      using XKGLDestroyRealFFTPlan() unless plan == NULL.
 *
 */
void KGLForwardRealFFT( //begin{proto}
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan **plan,   /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double complex *cdata,   /**< (out) returns the Fourier transform */
    const double   *rdata,   /**< (in) gives a real data */
    const size_t   size      /**< (in) gives the size of the real data */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,cdata != NULL,"cdata is NULL");
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
        plan1 = KGLCreateRealFFTPlan(status,KGL_FFT_FORWARD,size);
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
        if(!fftonce) *plan = plan1;
    }
    else {  /* !fftnew */
        plan1 = *plan;
        if(plan1->size != size) {
            KGLAddError(status,"plan->size != size\n");
            return;
        }

        r2c = plan1->r2c;
        if(r2c == NULL) {
            KGLAddRealFFTPlan(status,plan1);
            if(KGLCheckError(status)) return;
            r2c = plan1->r2c;
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

//    for(k = 0; k < size; k++) rdataP[k] = rdata[k];
    memcpy(rdataP,rdata,sizeof(double)*size);
#ifdef KGL_PTHREAD
    fftw_execute_r2c(r2c,rdataP,cdataP);
#else
    fftw_execute(r2c);
#endif
//    for(k = 0; k < size/2+1; k++) cdata[k] = cdataP[k];
    memcpy(cdata,cdataP,sizeof(fftw_complex)*(size/2+1));

#ifdef KGL_PTHREAD
    fftw_free(rdataP);
    fftw_free(cdataP);
#endif

    if(fftonce) KGLDestroyRealFFTPlan(plan1);

    return;
}

void KGLForwardRealFFTOnce( //begin{proto}
    KGLStatus      *status, /**< status */
    double complex *cdata,  /**< (out) returns the Fourier transform */
    const double   *rdata,  /**< (in) gives a real data */
    const size_t   size     /**< (in) gives the size of the real data */
    ) //end{proto}
{
    KGLForwardRealFFT(status,NULL,cdata,rdata,size);
    return;
}

/**
 * to calculate reverse FFT
 *
 * KGLRealFFTPlan **plan: see KGLForwardRealFFT()
 *
 */
void KGLReverseRealFFT( //begin{proto}
    KGLStatus      *status,  /**< status */
    KGLRealFFTPlan **plan,   /**< (in/out) KGLRealFFTPlan will be stored in *plan */
    double         *rdata,   /**< (out) returns the reverse Fourier transform */
    const double complex *cdata, /**< (in) gives a complex data */
    const size_t   size      /**< (in) gives the size of the real data */
    ) //end{proto}
{
    /* check on parameters */
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,cdata != NULL,"cdata is NULL");
    KGLAssert(status,rdata != NULL,"rdata is NULL");
    KGLAssert(status,size > 0,"size must be larger than 0");
    if(KGLCheckError(status)) return;
    /* end of check on parameters */

    bool fftonce = false;
    bool fftnew  = false;

    KGLRealFFTPlan *plan1  = NULL;
    double         *rdataP = NULL;
    fftw_complex   *cdataP = NULL;
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
        plan1 = KGLCreateRealFFTPlan(status,KGL_FFT_REVERSE,size);
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
        c2r     = plan1->c2r;
        if(!fftonce) *plan = plan1;
    }
    else { /* !fftwnew */
        plan1 = *plan;
        if(plan1->size != size) {
            KGLAddError(status,"plan->size != size\n");
            return;
        }

        c2r = plan1->c2r;
        if(c2r == NULL) {
            KGLAddRealFFTPlan(status,plan1);
            if(KGLCheckError(status)) return;
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
            KGLAddError(status,"fail to alloc rdata and/or cdata\n");
            return;
        }
#endif
    }

//    for(int k = 0; k < size/2+1; k++) cdataP[k] = cdata[k];
    memcpy(cdataP,cdata,sizeof(fftw_complex)*(size/2+1));
#ifdef KGL_PTHREAD
    fftw_execute_c2r(c2r,cdataP,rdataP);
#else
    fftw_execute(c2r);
#endif
//    for(int k = 0; k < size; k++) rdata[k] = rdataP[k];
    memcpy(rdata,rdataP,sizeof(double)*size);

#ifdef KGL_PTHREAD
    fftw_free(rdataP);
    fftw_free(cdataP);
#endif

    if(fftonce) KGLDestroyRealFFTPlan(plan1);

    return;
}

void KGLReverseRealFFTOnce( //begin{proto}
    KGLStatus      *status,  /**< status */
    double         *rdata,   /**< (out) returns the reverse Fourier transform */
    const double complex *cdata, /**< (in) gives a complex data */
    const size_t   size      /**< (in) gives the size of the real data */
    ) //end{proto}
{
    KGLReverseRealFFT(status,NULL,rdata,cdata,size);
    return;
}
