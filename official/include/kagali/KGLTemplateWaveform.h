#ifndef TEMPLATE_WAVEFORM_H
#define TEMPLATE_WAVEFORM_H

#include <complex.h>

KGL_BEGIN_DECLS

void KGLGenerateFreqPow( 
    int n,
    double fs,
    double power,
    double *freqpow
    );

void KGLGenerateTimePow( 
    int n,
    double fs,
    double power,
    double *timepow
    );

void KGLInspiralPhaseTD( 
    int iPN,
    double M,
    double eta,
    int n,
    double tmin,
    double tmax,
    double fs,
    double *phase
    );

void KGLInspiralPhaseFD( 
    int    iPN,
    double M,
    double eta,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double *phase
    );

void KGLGenerateInspiralTemplateFD( 
    int    iPN,
    int    npoint,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double M,
    double eta,
    double *freqpow,
    double complex *a1,
    double complex *a2
    );

void KGLGenerateInspiralSignalTD( 
    int iPN,
    int npoint,
    double tmin,
    double tmax,
    double fs,
    double M,
    double eta,
    double *timepow,
    double *a1,
    double *a2
    );

void KGLGenerateInspiralSignalFD( 
    int iPN,
    int npoint,
    double fmin,
    double fmax,
    double fs,
    double tc,
    double phic,
    double M,
    double eta,
    double *freqpow,
    double complex *a
    );

void KGLNormInspiral( 
    int npoint,
    double *Sn,
    double *freqpow,
    double fmin,
    double fmax,
    double fs,
    double M,
    double eta,
    double *norm
    );

void KGLMcEtaToTau03( 
    double Mc,
    double eta,
    double f0,
    double *tau0,
    double *tau3
    );

void KGLTau03ToMcEta( 
    double tau0,
    double tau3,
    double f0,
    double *Mc,
    double *eta
    );

void KGLTau03ToM12( 
    double tau0,
    double tau3,
    double f0,
    double *m1,
    double *m2
    );

void KGLTime2Freq3PN( 
    double M,
    double eta,
    double time,
    double *freq
    );


KGL_END_DECLS

#endif /* TEMPLATE_WAVEFORM_H */
