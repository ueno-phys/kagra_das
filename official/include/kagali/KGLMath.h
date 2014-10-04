#ifndef KGLMATH_H
#define KGLMATH_H

#include <math.h>

#define Icx _Complex_I

/**
 * \name Mathematical constants taken from LALConstants.h
 */
/*@{*/
#define KGL_E         2.7182818284590452353602874713526625  /**< e */
#define KGL_LOG2E     1.4426950408889634073599246810018922  /**< log_2 e */
#define KGL_LOG10E    0.4342944819032518276511289189166051  /**< log_10 e */
#define KGL_LN2       0.6931471805599453094172321214581766  /**< log_e 2 */
#define KGL_LN10      2.3025850929940456840179914546843642  /**< log_e 10 */
#define KGL_SQRT2     1.4142135623730950488016887242096981  /**< sqrt(2) */
#define KGL_SQRT1_2   0.7071067811865475244008443621048490  /**< 1/sqrt(2) */
#define KGL_GAMMA     0.5772156649015328606065120900824024  /**< gamma */
#define KGL_EXPGAMMA  1.7810724179901979852365041031071795  /**< exp(gamma) */
#define KGL_PI        3.1415926535897932384626433832795029  /**< pi */
#define KGL_TWOPI     6.2831853071795864769252867665590058  /**< 2*pi */
#define KGL_PI_2      1.5707963267948966192313216916397514  /**< pi/2 */
#define KGL_PI_4      0.7853981633974483096156608458198757  /**< pi/4 */
#define KGL_1_PI      0.3183098861837906715377675267450287  /**< 1/pi */
#define KGL_2_PI      0.6366197723675813430755350534900574  /**< 2/pi */
#define KGL_2_SQRTPI  1.1283791670955125738961589031215452  /**< 2/sqrt(pi) */
#define KGL_PI_180    1.7453292519943295769236907684886127e-2 /**< pi/180 */
#define KGL_180_PI    57.295779513082320876798154814105170 /**< 180/pi */
/*@}*/

KGL_BEGIN_DECLS

#ifdef __GNUC__
/* useful macros */
#define MAXd(a,b) (                                                     \
        {                                                               \
            double xkgl_a = (double)(a), xkgl_b = (double)(b);          \
            xkgl_a >= xkgl_b ? xkgl_a : xkgl_b;                         \
        }                                                               \
        )
#define MINd(a,b) (                                                     \
        {                                                               \
            double xkgl_a = (double)(a), xkgl_b = (double)(b);          \
            xkgl_a <= xkgl_b ? xkgl_a : xkgl_b;                         \
        }                                                               \
        )
#define MAXi(a,b) (                                                     \
        {                                                               \
            int xkgl_a = (int)(a), xkgl_b = (int)(b);                   \
            xkgl_a >= xkgl_b ? xkgl_a : xkgl_b;                         \
        }                                                               \
        )
#define MINi(a,b) (                                                     \
        {                                                               \
            int xkgl_a = (int)(a), xkgl_b = (int)(b);                   \
            xkgl_a <= xkgl_b ? xkgl_a : xkgl_b;                         \
        }                                                               \
        )
#define SQ(a) (                                 \
        {                                       \
            double xkgl_a = (a);                \
            xkgl_a*xkgl_a;                      \
        }                                       \
        )
#else
/* define functions for non-gcc */
static inline double MAXd(const double a, const double b) {
    return a >= b ? a : b;
}
static inline double MINd(const double a, const double b) {
    return a <= b ? a : b;
}
static inline int MAXi(const int a, const int b) {
    return a >= b ? a : b;
}
static inline int MINi(const int a, const int b) {
    return a <= b ? a : b;
}
static inline SQ(const double a) {
    return a*a;
}
#endif /* __GNUC__ */

KGL_END_DECLS

#endif /* KGLMATH_H */
