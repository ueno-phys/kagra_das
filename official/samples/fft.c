#include <kagali/KGLStdlib.h>
#include <kagali/RealFFT.h>
#include <kagali/BandPassFilter.h>

#include <sys/types.h>
#include <unistd.h>
#include <string.h>

typedef double complex cdouble;

#define NDATA_DEFAULT 4096

#define CALLOC(xx,nn,ss) do {                                   \
        if((xx = (ss *)calloc((nn),sizeof(ss))) == NULL) {      \
            KGLPrintError("fail to alloc \"%s\"",#xx);          \
            exit(EXIT_FAILURE);                                 \
        }                                                       \
    } while(0)

int main(int argc, char *argv[]) {
    KGLStatus *status = KGLCreateStatus();

    size_t ndata = NDATA_DEFAULT;       /* the size of data */

    if(argc > 2) {
        fprintf(stderr,"usage: %s [size[k|m]]\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    else if(argc == 2) {
        char *str = argv[1];
        char *endptr = NULL;
        long int nd = strtol(argv[1],&endptr,0);
        if(endptr == str) {
            fprintf(stderr,"No digits are given.\n");
            exit(EXIT_FAILURE);
        }
        if(nd <= 0) {
            if(strlen(endptr) != 1) {
                fprintf(stderr,"Invalid number is given.\n");
                exit(EXIT_FAILURE);
            }
        }
        if(*endptr != '\0') {
            if(strlen(endptr) != 1) {
                fprintf(stderr,"Invalid number is given.\n");
                exit(EXIT_FAILURE);
            }
            switch(*endptr) {
            case 'k':
            case 'K':
                ndata = nd*1024;
                break;
            case 'm':
            case 'M':
                ndata = nd*1024*1024;
                break;
            default:
                fprintf(stderr,"Invalid number is given.\n");
                exit(EXIT_FAILURE);
                break;
            }
        }
        else {
            ndata = nd;
        }
    }
    printf("data size = %lu\n",ndata);

    size_t nfreq = (ndata/2+1);

    double t0   = 0.0;                /* the start time (in sec) */
    double Tlen = 2.0;                /* the length of signal (in sec) */
    double dt   = Tlen/(double)ndata; /* sampling interval */
    double df   = 1.0/Tlen;

    FILE *fout = NULL;
    char *fname = NULL;

    KGLRealFFTPlan *plan = NULL;

    double *tsec = NULL, *sig = NULL, *freq = NULL;
    cdouble *sigf = NULL;

    CALLOC(tsec,ndata,double);
    CALLOC(sig,ndata,double);
    CALLOC(freq,nfreq,double);
    CALLOC(sigf,nfreq,cdouble);

    for(int k = 0; k < ndata; k++) tsec[k] = t0 + dt*k;
    for(int k = 0; k < nfreq; k++) freq[k] = df*k;

    /*     freq (Hz) , amplitude, phase */
    double f1 = 100.0,  a1 = 1.0, d1 = 0.0;
    double f2 = 200.25, a2 = 1.5, d2 = 0.15*KGL_PI;
    double f3 = 300.0,  a3 = 2.0, d3 = KGL_PI/3.0;
    for(int k = 0; k < ndata; k++) {
        double ti = tsec[k];
        sig[k] = a1*cos(2*KGL_PI*f1*ti + d1)
               + a2*cos(2*KGL_PI*f2*ti + d2)
               + a3*cos(2*KGL_PI*f3*ti + d3);
    }

    KGLForwardRealFFT(status,&plan,sigf,sig,ndata);
    KGLAbortIfError(status);

    fname = "signal.data";
    if((fout = fopen(fname,"w")) == NULL) {
        perror(fname);
        abort();
    }
    for(int k = 0; k < ndata; k++) {
        fprintf(fout,"%18.10e %18.10e\n",tsec[k],sig[k]);
    }
    fclose(fout);

    fname = "spect.data";
    if((fout = fopen(fname,"w")) == NULL) {
        perror(fname);
        abort();
    }
    for(int k = 0; k < nfreq; k++) {
        double ps = 2*cabs(sigf[k])/(double)ndata;
        fprintf(fout,"%18.10e %18.10e\n",freq[k],ps);
    }
    fclose(fout);

    KGLBandPassFilter(status,&plan,sig,ndata,dt,200.0,500.0);
    KGLAbortIfError(status);

    fname = "signal2.data";
    if((fout = fopen(fname,"w")) == NULL) {
        perror(fname);
        abort();
    }
    for(int k = 0; k < ndata; k++) {
        fprintf(fout,"%18.10e %18.10e\n",tsec[k],sig[k]);
    }
    fclose(fout);

    KGLForwardRealFFT(status,&plan,sigf,sig,ndata);
    KGLAbortIfError(status);

    fname = "spect2.data";
    if((fout = fopen(fname,"w")) == NULL) {
        perror(fname);
        abort();
    }
    for(int k = 0; k < nfreq; k++) {
        double ps = 2*cabs(sigf[k])/(double)ndata;
        fprintf(fout,"%18.10e %18.10e\n",freq[k],ps);
    }
    fclose(fout);

    KGLDestroyRealFFTPlan(plan);
//    fftw_cleanup();

    free(tsec);
    free(freq);
    free(sig);
    free(sigf);
    KGLDestroyStatus(status);

    return 0;
}
