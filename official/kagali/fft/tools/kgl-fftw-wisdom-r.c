#include <kagali/KGLStdlib.h>
#include "../src/configRealFFT.h"
#include <kagali/RealFFT.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <fcntl.h>

#define DEFAULT_OUTPUT "wisdom-r"
#define DEFAULT_FROM 2
#define DEFAULT_TO   20

char *cmdname = NULL;
char *tmpoutput = NULL;

void usage();
void sigterm(int sig);
char *maketmpfile(char *str);
int mkwisdom(int p);

void usage() {
    fprintf(stderr,"usage: %s [options] [num1[:num2]]",cmdname);
    fprintf(stderr,"  Create wisdom for specified sizes\n\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"   -o FILE, --output=FILE: output to FILE (default: %s)\n",DEFAULT_OUTPUT);
    fprintf(stderr,"   -w FILE, --wisdom=FILE: read wisdom from FILE\n");
    fprintf(stderr,"   -u, --update: read wisdom from output (or wisdom) and update it.\n");
    fprintf(stderr,"                 Ignore this option if both output and wisdom are specified.\n");
    fprintf(stderr,"   -h, --help: display this message\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"num1[:num2]: create wisdom for 2^{num1} [or from 2^{num1} to 2^{num2}]\n");
    fprintf(stderr,"           (default: %d:%d)\n",DEFAULT_FROM,DEFAULT_TO);

    return;

}

void sigterm(int sig) {
    unlink(tmpoutput);
    fprintf(stderr,"\n%s: stop by signal: %d\n",cmdname,sig);
    exit(sig);
}

char *maketmpfile(char *str) {
    char *na = 
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    char *name;
    size_t name_length;
    char *c0, *c1, *c2, *c3, *c4;
    int fd;

    name_length = strlen(str) + 10;
    name = (char *)malloc(name_length);
    if(name == NULL) return NULL;

    for(c0 = na; *c0 != '\0'; c0++) {
        for(c1 = na; *c1 != '\0'; c1++) {
            for(c2 = na; *c2 != '\0'; c2++) {
                for(c3 = na; *c3 != '\0'; c3++) {
                    for(c4 = na; *c4 != '\0'; c4++) {
                        snprintf(name,name_length,"%s.%c%c%c%c%c",
                                 str,*c0,*c1,*c2,*c3,*c4);
                        fd = open(name,O_RDONLY);
                        if(fd == -1) {
                            if(errno == ENOENT) return name;
                        }
                        else {
                            close(fd);
                        }
                    }
                }
            }
        }
    }

    return NULL;
}

int mkwisdom(int p) {
    uint32_t flag = KGL_FFTWPLAN;
    int n = (1<<p);
    int retval = 0;
    double *rdata = NULL;
    fftw_complex *cdata = NULL;
    fftw_plan r2c = NULL, c2r = NULL;

    fprintf(stderr,"n = %d (2^%d):",n,p);
    fflush(stderr);

    rdata = (double *)fftw_malloc(sizeof(double)*n);
    if(rdata == NULL) {
        fprintf(stderr,"fail to malloc(rdata)\n");
        retval = -1;
        goto ret;
    }
    cdata = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(n/2+1));
    if(cdata == NULL) {
        fprintf(stderr,"fail to malloc(cdata)\n");
        retval = -1;
        goto ret;
    }

    fprintf(stderr," r2c");
    fflush(stderr);
    r2c = fftw_plan_dft_r2c_1d(n,rdata,cdata,flag);
    if(r2c == NULL) {
        retval = -1;
        goto ret;
    }

    fprintf(stderr," c2r");
    fflush(stderr);
    c2r = fftw_plan_dft_c2r_1d(n,cdata,rdata,flag);
    if(c2r == NULL) {
        retval = -1;
        goto ret;
    }

ret:
    fprintf(stderr,"\n");

    fftw_destroy_plan(r2c);
    fftw_destroy_plan(c2r);
    fftw_free(rdata);
    fftw_free(cdata);

    return retval;
}

int main(int argc, char *argv[]) {
    struct option long_opts[] = {
        {"output",1,0,'o'},
        {"wisdom",1,0,'w'},
        {"update",0,0,'u'},
        {"help"  ,0,0,'h'},
        {0,0,0,0}
    };

    bool update = false;
    size_t rsize, wsize;
    int from = 0, to = 0;
    char *output = NULL;
    char *wisdom = NULL;
    struct stat st;
    char buf[BUFSIZ];

    FILE *wfile, *ofile;

    fprintf(stderr,"Wwisdom will be created for %s\n",fftw_version);

    cmdname = strrchr(argv[0],'/');
    if(cmdname == NULL) {
        cmdname = argv[0];
    }
    else {
        cmdname++;
        if(*cmdname == '\0') cmdname = argv[0];
    }

    while(1) {
        int option_index = 0;
        int c = getopt_long(argc,argv,"o:w:uh",long_opts,&option_index);
        if(c == -1) break;

        switch(c) {
        case 'o':
            output = optarg;
            break;
        case 'w':
            wisdom = optarg;
            break;
        case 'u':
            update = true;
            break;
        case 'h':
            usage();
            exit(0);
            break;
        case '?':
            fprintf(stderr,"unknown option %c\n",optopt);
            exit(EXIT_FAILURE);
            break;
        default:
            fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
            exit(EXIT_FAILURE);
            break;
        }
    }

    switch(argc-optind) {
    case 0:
        from = DEFAULT_FROM;
        to = DEFAULT_TO;
        break;
    case 1:
        do {
            char *argv1, *endp;
            bool error = false;
            argv1 = argv[optind];
            from = strtod(argv1,&endp);
            to = from;
            if(!(error = (argv1 == endp))) {
                    if(*endp == '\0') {
                        to = from;
                    }
                    else if(*endp == ':') {
                        argv1 = endp + 1;
                        to = strtod(argv1,&endp);
                        error = (argv1 == endp);
                    }
            }
            if(error) {
                fprintf(stderr,"invalid argment num1 and/or num2\n\n");
                usage();
                exit(EXIT_FAILURE);
            }
        } while(0);
        break;
    default:
        usage();
        exit(EXIT_FAILURE);
    }

    if(update) {
        if(output == NULL) {
            if(wisdom != NULL) {
                output = wisdom;
            }
        }
        if(wisdom == NULL) {
            if(output != NULL) {
                wisdom = output;
            }
            else {
                wisdom = DEFAULT_OUTPUT;
            }
        }
    }

    if(output == NULL) output = DEFAULT_OUTPUT;

    if(wisdom != NULL) {
        wfile = fopen(wisdom,"r");
        if(wfile == NULL) {
            if(errno != ENOENT) {
                perror(wisdom);
                abort();
            }
        }
        else {
            if(stat(wisdom,&st) != 0) {
                perror(wisdom);
                abort();
            }
            if(st.st_size > 0) {
                int ret = fftw_import_wisdom_from_file(wfile);
                if(ret == 0) {
                    fprintf(stderr,"fail to import wisdom from %s\n",wisdom);
                    abort();
                }
            }
            fclose(wfile);
        }
    }

    tmpoutput = maketmpfile(output);
    if(tmpoutput == NULL) {
        fprintf(stderr,"fail to tmpoutput name\n");
        abort();
    }

    wfile = fopen(tmpoutput,"w");
    if(wfile == NULL) {
        perror(output);
        abort();
    }
    signal(SIGINT,sigterm);
    signal(SIGTERM,sigterm);

    for(int p = from; p <= to; p++) {
        if(mkwisdom(p) < 0) abort();
    }

    fftw_export_wisdom_to_file(wfile);
    if(fclose(wfile) != 0) {
        perror(tmpoutput);
        abort();
    }

    signal(SIGINT,SIG_DFL);
    signal(SIGTERM,SIG_DFL);

    wfile = fopen(tmpoutput,"r");
    if(wfile == NULL) {
        perror(tmpoutput);
        exit(EXIT_FAILURE);
    }

    ofile = fopen(output,"w");
    if(ofile == NULL) {
        perror(output);
        exit(EXIT_FAILURE);
    }

    while((rsize = fread(buf,1,BUFSIZ,wfile)) > 0) {
        wsize = fwrite(buf,1,rsize,ofile);
        if(wsize <= 0) {
            fprintf(stderr,"Error occured in writing output to %s\n",output);
            abort();
        }
    }
    if(ferror(wfile)) {
        fprintf(stderr,"Error occured in reading file %s\n",tmpoutput);
        abort();
    }

    fclose(ofile);
    fclose(wfile);

    unlink(tmpoutput);

    return 0;
}
