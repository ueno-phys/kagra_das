#ifndef NOISE_PSD_H
#define NOISE_PSD_H

#include <kagali/KGLStdlib.h>

KGL_BEGIN_DECLS

void KGLNoiseGEO( 
    int    n,
    double *Sn,
    double fs
    );

void KGLNoiseLIGOI( 
    int    n,
    double *Sn,
    double fs
    );

void KGLNoiseLIGOA( 
    int    n,
    double *Sn,
    double fs
    );

void KGLNoiseVIRGO( 
    int    n,
    double *Sn,
    double fs
    );

void KGLReadNoiseSpectrum( 
    int    iDet,
    int    n,
    double *Sn,
    double fs
    );


KGL_END_DECLS

#endif /* NOISE_PSD_H */
