#include <kagali/KGLStdlib.h>

void KGLtest1(KGLStatus *status,size_t size);

void KGLtest2(KGLStatus *status,
           double *outarray, double *inarray, size_t size);

void KGLtest1(KGLStatus *status, size_t size) {
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,size > 0,"size must be larger than 0");
    if(KGLCheckError(status)) return;

    double *vout = NULL;
    double *vin  = NULL;
    vout = (double *)calloc(size,sizeof(double));
    vin  = (double *)calloc(size,sizeof(double));
    if((vout == NULL) || (vin == NULL)) {
        KGLAddError(status,"fail to alloc vout and/or vin\n");
        goto free_return;
    }

    KGLtest2(status,vout,vin,size);
    if(KGLCheckError(status)) {
        KGLAddMessage(status,"");
        goto free_return;
    }

free_return:
    free(vout);
    free(vin);

    return;
}

void KGLtest2(KGLStatus *status,
           double *outarray, double *inarray, size_t size)
{
    if(status == NULL) {
        KGLPrintError("NULL status pointer passed");
        abort();
    }

    KGLAssert(status,outarray != NULL,"outarray is NULL");
    KGLAssert(status,inarray  != NULL,"inarray is NULL");
    KGLAssert(status,size > 0,"size must be larger than 0");
//    KGLAssert(status,size > 0,NULL);
//    KGLAssert(status,size > 0,"");
    if(KGLCheckError(status)) return;

    for(int k = 0; k < size; k++) outarray[k] = inarray[k];

    return;
}

int main(int argc, char *argv[]) {
    KGLStatus *status = KGLCreateStatus();

    size_t size = 10;
//    double *in = NULL;
//    double *out = NULL;
//    double in[size], out[size];

//    KGLtest2(status,out,in,size);
    KGLtest1(status,size);
    KGLAbortIfError(status);

    KGLDestroyStatus(status);

    return 0;
}
