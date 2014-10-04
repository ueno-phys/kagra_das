/**
 * \author K. Oohara
 *
 * \brief Includes the standard KAGALI header files.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <kagali/KGLStdlib.h>
 * \endcode
 *
 * This header is the overall header for the \c std
 * package.
 *
 * This header also includes function prototype headers for certain standard modules used
 * by many KAGALI routines:
 *
 * \code
 * #include <stdio.h>
 * #include <stdlib.h>
 * #include <stdarg.h>
 * #include <stdbool.h>
 * #include <inttypes.h>
 * #include <kagali/KGLError.h>
 * #include <kagali/KGLMath.h>
 * \endcode
 *
 */

#ifndef KGLSTDLIB_H
#define KGLSTDLIB_H

#ifdef  __cplusplus
# define KGL_BEGIN_DECLS  extern "C" {
# define KGL_END_DECLS    }
#else
# define KGL_BEGIN_DECLS
# define KGL_END_DECLS
#endif

KGL_BEGIN_DECLS
volatile char *GITVER_KGL;
KGL_END_DECLS

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <inttypes.h>
#include <errno.h>

#include <kagali/KGLError.h>
#include <kagali/KGLMath.h>

/* #undef KGL_PTHREAD_LOCK */

#endif /* KGLSTDLIB_H */
