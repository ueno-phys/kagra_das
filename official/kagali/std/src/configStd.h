#ifndef CONFIGSTD_H
#define CONFIGSTD_H

#ifdef KGL_PTHREAD
# define KGL_STD_LOCK   XKGLPthreadLock()
# define KGL_STD_UNLOKC XKGLPthreadUnlock()
#else
# define KGL_STD_LOCK   /* nothing */
# define KGL_STD_UNLOCK /* nothing */
#endif

#endif /* CONFIGSTD_H */
