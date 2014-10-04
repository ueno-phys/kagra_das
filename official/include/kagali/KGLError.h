#ifndef KGLERROR_H
#define KGLERROR_H

KGL_BEGIN_DECLS

typedef struct tagKGLStatus {
    int statusCode;
    size_t messageSize;
    char *message;
} KGLStatus;

#define KGL_EL_NONE    (0)
#define KGL_EL_INFO    (1 << 0)
#define KGL_EL_WARNING (1 << 1)
#define KGL_EL_ERROR   (1 << 2)

#ifdef KGL_GITID
# define KGL_FILENAME __FILE__ "(" KGL_GITID ")"
#else
# define KGL_FILENAME __FILE__
#endif

/**
   KGL_StatusVal() extract the statusCode.
   This can be used as the lvalue.
 */
#define KGL_StatusVal(status) (status)->statusCode

/**
   KGLCheckStatus() and KGLCheckError() check whether error occured.
 */
#define KGLCheckStatus(status,level) (KGL_StatusVal(status) >= level)
#define KGLCheckError(status) KGLCheckStatus(status,KGL_EL_ERROR)

/**
   KGLAddToStatus() adds a message to status.
   The value of code must be one of KGL_EL_NONE, KGL_EL_INFO, KGL_EL_WARNING
   and KGL_EL_ERROR.
 */
#define KGLAddToStatus(status,code,...)                               \
    XKGLAddToStatus(status,code,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)
#define KGLAddError(status,...)        \
    XKGLAddToStatus(status,KGL_EL_ERROR,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)
#define KGLAddWarning(status,...)        \
    XKGLAddToStatus(status,KGL_EL_WARNING,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)
#define KGLAddInfo(status,...)        \
    XKGLAddToStatus(status,KGL_EL_INFO,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)
#define KGLAddMessage(status,...)        \
    XKGLAddToStatus(status,KGL_EL_NONE,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)

/**
   initialization of status
 */
#define KGLCreateStatus()                               \
    XKGLCreateStatus(__func__,KGL_FILENAME,__LINE__)

#define KGLInitStatus(status)                                           \
    XKGLInitStatus(status,__func__,KGL_FILENAME,__LINE__)

#define KGLResetStatus(status)                                           \
    XKGLResetStatus(status,__func__,KGL_FILENAME,__LINE__)

#define KGLPrintStatusMessage(status)                                   \
    XKGLPrintStatusMessage(status,__func__,KGL_FILENAME,__LINE__)

#define KGLPrintMsg(header,...)                                         \
    XKGLPrintMsg(header,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)

#define KGLPrintError(...)                                              \
    XKGLPrintMsg("Error",__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)

#define KGLPrintWarning(...)                                            \
    XKGLPrintMsg("Warning",__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)

#define KGLPrintInfo(...)                                       \
    XKGLPrintMsg(NULL,__func__,KGL_FILENAME,__LINE__,__VA_ARGS__)

#if defined KGL_NDEBUG || defined NDEBUG
# define KGLAssert(status,expr,msg) /* nothing */
#else
# define KGLAssert(status,expr,msg) do {                                \
        if(expr) {                                                      \
        }                                                               \
        else {                                                          \
            XKGLAssert(status,                                          \
                       msg,                                             \
                       "Assertion: `"__STRING(expr)"' failed",          \
                       __func__,KGL_FILENAME,__LINE__);                 \
        }                                                               \
    } while(0)
#endif

#define KGLAbort(status) XKGLAbort(status,__func__,KGL_FILENAME,__LINE__)
#define KGLAbortIfError(status) if(KGLCheckError(status)) KGLAbort(status)

void XKGLAddToStatus( 
    KGLStatus *status,
    const int code,
    const char *func, const char *file, int line,
    const char *fmt,
    ...
    );

KGLStatus *XKGLCreateStatus( 
    const char *func, const char *file, int line
    );

void XKGLInitStatus( 
    KGLStatus *status,
    const char *func, const char *file, int line
    );

void XKGLResetStatus( 
    KGLStatus *status,
    const char *func, const char *file, int line
    );

void KGLDestroyStatus(KGLStatus *status);

void XKGLVPrintMsg( 
    const char *header,
    const char *func, const char *file, int line,
    const char *fmt,
    va_list ap
    );

void XKGLPrintMsg( 
    const char *header,
    const char *func, const char *file, int line,
    const char *fmt,...
    );

void XKGLAssert( 
    KGLStatus *status,
    const char *msg,
    const char *expr,
    const char *func, const char *file, int line
    );

void XKGLPrintStatusMessage( 
    const KGLStatus *status,
    const char *func, const char *file, int line
    );

void XKGLAbort( 
    const KGLStatus *status,
    const char *func, const char *file, int line
    );


KGL_END_DECLS

#endif /* KGLERROR_H */
