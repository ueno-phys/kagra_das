#ifndef LS_FUNC_H
#define LS_FUNC_H

KGL_BEGIN_DECLS

void KGLCostLS( 
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double *F
    );

void KGLDerivCostLS( 
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double *dFdA,
    double *dFdf,
    double *dFdp
    );

void KGLHessian3D( 
    double fs,
    int nframe,
    double *frame,
    double A,
    double f,
    double phi,
    double **H
    );

void KGLInvMatrix( 
    int n,
    double **H
    );


KGL_END_DECLS

#endif /* LS_FUNC_H */
