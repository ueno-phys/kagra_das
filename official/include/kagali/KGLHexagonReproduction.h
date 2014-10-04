#ifndef HEX_REPRODUCTION_H
#define HEX_REPRODUCTION_H

KGL_BEGIN_DECLS

void KGLReprodHexagonGrids( 
    int    imode,    /**< (1 or 2) choose one of the two eigen vectors */
    double mm,       /**< minimal match */
    double lambda1,  /**< eigen value 1 of the metric */
    double lambda2,  /**< eigen value 2 of the metric */
    double theta,    /**< inclination angle of ellipse */
    double tau0,     /**< mass parameter which is used to define the metric */
    double tau3,     /**< mass parameter which is used to define the metric */
    double *tau0d,   /**< (out) mass parameter grids which are newly generated */
    double *tau3d    /**< (out) mass parameter grids which are newly generated */
    );

void KGLReprodHexagonVertices( 
    int    imode,    /**< (1 or 2) choose one of the two eigen vectors */
    double mm,       /**< minimal match */
    double lambda1,  /**< eigen value 1 of the metric */
    double lambda2,  /**< eigen value 2 of the metric */
    double theta,    /**< inclination angle of ellipse */
    double tau0,     /**< mass parameter which is used to define the metric */
    double tau3,     /**< mass parameter which is used to define the metric */
    double *tau0v,   /**< (out) Hexagon vertices around a given (tau0,tau3) */
    double *tau3v    /**< (out) Hexagon vertices around a given (tau0,tau3) */
    );


KGL_END_DECLS

#endif /* HEX_REPRODUCTION_H */
