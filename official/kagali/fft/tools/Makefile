CC := gcc
GSLLIBS := -lgsl -lgslcblas
FFTLIBS := -lfftw3
INCDIR :=
pkgdir := ../../..
pkgincdir := $(pkgdir)/include
pkglibdir := $(pkgdir)/lib
scriptdir := $(pkgdir)/scripts

PROGS := \
	kgl-fftw-wisdom-r

OBJS := ${PROGS:%=%.o}

CFLAGS  := -ggdb -std=c99 -Wall -O2 $(INCDIR) -I$(pkgincdir)
LDFLSGS := -ggdb -std=c99 -Wall -O2

$(PROGS):

kgl-fftw-wisdom-r: kgl-fftw-wisdom-r.o 
	$(CC) -o $@ $^ $(LDFLAGS) $(FFTLIBS)

kgl-fftw-wisdom-r.o: kgl-fftw-wisdom-r.c ../src/configRealFFT.h
	$(CC) -c $< $(CFLAGS)

test:
	$(CC) $(CFLAGS) -E kgl-fftw-wisdom-r.c

clean:
	rm -f $(PROGS) $(OBJS)
