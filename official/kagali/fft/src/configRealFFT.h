#ifndef CONFIGREALFFT_H
#define CONFIGREALFFT_H

#define KGL_FFTWPLAN   (FFTW_PATIENT|FFTW_DESTROY_INPUT)
#define KGL_FFTWPLAN_FS (FFTW_ESTIMATE|FFTW_DESTROY_INPUT)

#ifdef KGL_PTHREAD
# define KGL_FFTW_LOCK   XKGLPthreadLock()
# define KGL_FFTW_UNLOKC XKGLPthreadUnlock()
#else
# define KGL_FFTW_LOCK   /* nothing */
# define KGL_FFTW_UNLOCK /* nothing */
#endif

#endif /* CONFIGREALFFT_H */
