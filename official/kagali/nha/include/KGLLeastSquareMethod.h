#ifndef LS_METHOD_H
#define LS_METHOD_H

KGL_BEGIN_DECLS

void KGLSteepestDescent( 
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double mu,
    double A,
    double *f,
    double *phi,
    double *Fcost
    );

void KGLNewton3D( 
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double mu,
    double *A,
    double *f,
    double *phi
    );

void KGLAmpConv( 
    double fs,
    int nframe,
    double *frame,
    int nitr,
    double nu,
    double f,
    double phi,
    double *A,
    double *Fcost
    );


KGL_END_DECLS

#endif /* LS_METHOD_H */
