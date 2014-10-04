#include "KGLVersion.h"
#include <kagali/KGLStdlib.h>
#include "configStd.h"
#include <string.h>

void XKGLAddToStatus( //begin{proto}
    KGLStatus *status,
    const int code,
    const char *func, const char *file, int line,
    const char *fmt,
    ...
    ) //end{proto}
{
    va_list ap, apb;
    int msgLen, headerLen;
    size_t oldSize, newSize;

    char msgtest[10];
    char header[100];

    char *oldmessage = NULL;
    char *newmessage = NULL;
    char *msgptr = NULL;

    const bool funcNotNull = ((func != NULL) && (*func != '\0'));
    const char *funcfmt = " < in function %s (line %d of %s)\n";

    if(status == NULL) {
        XKGLPrintMsg("Error",func,file,line,"NULL status pointer passed");
        abort();
    }

    KGL_STD_LOCK;

    if(code != 0) status->statusCode |= code;

    oldmessage = status->message;
    oldSize    = status->messageSize;

    header[0] = '\0';
    switch (code){
    case KGL_EL_NONE:
        break;
    case KGL_EL_INFO:
        strncat(header,"KGL Info: ",sizeof(header)-1);
        break;
    case KGL_EL_WARNING:
        strncat(header,"KGL Waring: ",sizeof(header)-1);
        break;
    case KGL_EL_ERROR:
        strncat(header,"KGL Error: ",sizeof(header)-1);
        break;
    default:
        XKGLPrintMsg("Error",func,file,line,"invalid number is given for code");
        KGL_STD_UNLOCK;
        XKGLAbort(status,NULL,NULL,0);
        break;
    }
    headerLen = strlen(header);

    va_start(ap,fmt);
    va_copy(apb,ap);
    msgLen = vsnprintf(msgtest,sizeof(msgtest),fmt,ap);
    va_copy(ap,apb);

    if(funcNotNull) {
        msgLen += snprintf(msgtest,sizeof(msgtest),funcfmt,
                           func,line,file);
    }
    else {
        msgLen += snprintf(msgtest,sizeof(msgtest),"\n");
    }

    newSize = oldSize + headerLen + msgLen;
    newmessage = (char *)realloc(oldmessage,newSize+1);
    if(newmessage == NULL) {
        XKGLVPrintMsg(NULL,func,file,line,fmt,ap);
    }            
    else {
        if(headerLen > 0) strcat(newmessage,header);
        msgptr = newmessage + oldSize + headerLen;
        msgLen = vsprintf(msgptr,fmt,ap);
        msgptr += msgLen;
        if(funcNotNull) {
            sprintf(msgptr,funcfmt,func,line,file);
        }
        else {
            sprintf(msgptr,"\n");
        }
        status->message = newmessage;
        status->messageSize = strlen(status->message);
    }

    va_end(ap);
    va_end(apb);

    KGL_STD_UNLOCK;

    return;
}

KGLStatus *XKGLCreateStatus( //begin{proto}
    const char *func, const char *file, int line
    ) //end{proto}
{
    KGLStatus *status = (KGLStatus *)calloc(1,sizeof(KGLStatus));
    if(status == NULL) {
        XKGLPrintMsg("Error",func,file,line,"Fail to alloc status");
        abort();
    }
    XKGLInitStatus(status,func,file,line);

    return status;
}

void XKGLInitStatus( //begin{proto}
    KGLStatus *status,
    const char *func, const char *file, int line
    ) //end{proto}
{
    if(status == NULL) {
        XKGLPrintMsg("Error",func,file,line,"NULL status pointer passed");
        abort();
    }

    KGL_STD_LOCK;

    status->statusCode = KGL_EL_NONE;
    status->messageSize = 0;
    status->message = (char *)malloc(1);
    if(status->message == NULL) {
        XKGLPrintMsg("Error",func,file,line,
                     "fail to allocate memory for status messages");
        KGL_STD_UNLOCK;
        abort();
    }
    status->message[0] = '\0';

    KGL_STD_UNLOCK;

    return;
}

void XKGLResetStatus( //begin{proto}
    KGLStatus *status,
    const char *func, const char *file, int line
    ) //end{proto}
{
    if(status == NULL) {
        XKGLPrintMsg("Error",func,file,line,"NULL status pointer passed");
        abort();
    }

    KGL_STD_LOCK;

    status->statusCode = KGL_EL_NONE;
    status->messageSize = 0;
    free(status->message);
    status->message = (char *)malloc(1);
    if(status->message == NULL) {
        XKGLPrintMsg("Error",func,file,line,
                     "fail to allocate memory for status messages");
        KGL_STD_UNLOCK;
        abort();
    }
    status->message[0] = '\0';

    KGL_STD_UNLOCK;

    return;
}

void KGLDestroyStatus(KGLStatus *status) //prototype
{
    if(status == NULL) return;

    KGL_STD_LOCK;

    free(status->message);
    free(status);

    KGL_STD_UNLOCK;

    return;
}

/* Prints an error, warning or information messages */
void XKGLVPrintMsg( //begin{proto}
    const char *header,
    const char *func, const char *file, int line,
    const char *fmt,
    va_list ap
    ) //end{proto}
{
    if(header != NULL && *header != '\0') {
        fprintf(stderr,"KGL %s: ",header);
    }

    vfprintf(stderr,fmt,ap);

    if(func != NULL && *func != '\0') {
        fprintf(stderr," in %s (line %d of %s)\n",
                func,line,file);
    }
    else {
        fprintf(stderr,"\n");
    }
}

void XKGLPrintMsg( //begin{proto}
    const char *header,
    const char *func, const char *file, int line,
    const char *fmt,...
    ) //end{proto}
{
    va_list ap;
    va_start(ap,fmt);

    XKGLVPrintMsg(header,func,file,line,fmt,ap);

    va_end(ap);

    return;
}

void XKGLAssert( //begin{proto}
    KGLStatus *status,
    const char *msg,
    const char *expr,
    const char *func, const char *file, int line
    ) //end{proto}
{
    if((msg == NULL) || (*msg == '\0')) {
        XKGLAddToStatus(status,KGL_EL_ERROR,func,file,line,"%s\n",expr);
    }
    else {
        XKGLAddToStatus(status,KGL_EL_ERROR,func,file,line,"Assertion: %s\n",
                          msg);
    }
    return;
}

void XKGLPrintStatusMessage( //begin{proto}
    const KGLStatus *status,
    const char *func, const char *file, int line
    ) //end{proto}
{
    if(status == NULL) {
        XKGLPrintMsg("Warning",func,file,line,"NULL status pointer passed");
    }
    else {
        fprintf(stderr,"%s",status->message);
    }

    return;
}

void XKGLAbort( //begin{proto}
    const KGLStatus *status,
    const char *func, const char *file, int line
    ) //end{proto}
{
    XKGLPrintStatusMessage(status,func,file,line);
    fprintf(stderr,"KGL Abort");
    if(func != NULL && *func != '\0') {
        fprintf(stderr," in function %s (line %d of %s)",func,line,file);
    }
    fprintf(stderr,"\n");

    abort();
}
