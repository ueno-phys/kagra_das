#ifndef TEMPLATE_METRIC_H
#define TEMPLATE_METRIC_H

#include <complex.h>

KGL_BEGIN_DECLS

void KGLNoiseMoment( 
   int    iq,
   int    ir,
   double *Sn,
   int    n,
   double fmin,
   double fmax,
   double fs,
   double *val);

void KGLNoiseMoment7( 
    double *Sn,
    int    n,
    double fmin,
    double fmax,
    double fs,
    double **vr
    );

void KGLPhaseGradExpansionCoeff( 
    double tau0,
    double tau3,
    double fmin,
    double **psi
    );

void KGLTemplateMetric( 
    double *Sn,
    int n,
    double tau0,
    double tau3,
    double fmin,
    double fmax,
    double fs,
    double **vr,
    double **gml
    );

void KGLTemplateMetricEigen( 
    double **gml,
    double *theta,
    double *lambda1,
    double *lambda2
    );

void KGLTemplateMetricDistance( 
    double **gml,
    double dtau0,
    double dtau3,
    double *ds2
    );


KGL_END_DECLS

#endif /* TEMPLATE_METRIC_H */
