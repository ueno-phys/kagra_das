GSLLIBS := -lgsl -lgslcblas
FFTLIBS := -lfftw3

KGLDIR  := /Users/kagra/official
KGLINCS := -I$(KGLDIR)/include -I/opt/local/include
KGLLIBS := -L$(KGLDIR)/lib -lkagali

CFLAGS := -ggdb -O2 -std=c99  -Wall $(KGLINCS)
LDFLAGS := -ggdb -O2 -std=c99 -Wall $(KGLLIBS)

SRCS := \
	sample1.c \
        fft.c \
	KGLGenerateHexagonalTemplateBank.c \
	nha.c

OBJS := ${SRCS:%.c=%}

all: $(OBJS)

sample1: sample1.o $(KGLDIR)/lib/libkagali.a
	$(CC) -I/opt/local/include -L/opt/local/lib -o $@ $^ $(LDFLAGS) $(GSLLIBS) $(FFTLIBS)

sample1-n: sample1-n.o $(KGLDIR)/lib/libkagali.a
	$(CC) -I/opt/local/include -L/opt/local/lib -o $@ $^ $(LDFLAGS) $(GSLLIBS) $(FFTLIBS)

sample1-n.o: sample1.c
	$(CC) -I/opt/local/include -c -o $@ $^ $(CFLAGS) -DKGL_NDEBUG

fft: fft.o $(KGLDIR)/lib/libkagali.a
	$(CC) -L/opt/local/lib -o $@ $^ $(LDFLAGS) $(GSLLIBS) $(FFTLIBS)

KGLGenerateHexagonalTemplateBank: KGLGenerateHexagonalTemplateBank.o $(KGLDIR)/lib/libkagali.a
	$(CC) -I/opt/local/include -L/opt/local/lib -o $@ $^ $(LDFLAGS) $(GSLLIBS) $(FFTLIBS)

nha: nha.o $(KGLDIR)/lib/libkagali.a
	$(CC) -I/opt/local/include -L/opt/local/lib -o $@ $^ $(LDFLAGS) $(GSLLIBS) $(FFTLIBS)

clean:
	rm -rf $(OBJS) *.o *.dSYM
