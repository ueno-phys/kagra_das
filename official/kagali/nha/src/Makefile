CC := gcc
pkgdir := ../../..
pkgincdir := $(pkgdir)/include
pkglibdir := $(pkgdir)/lib
scriptdir := $(pkgdir)/scripts
incdir := ../include

libdir := ../lib

SRCS := \
	KGLLeastSquareFunc.c \
	KGLLeastSquareMethod.c \
	KGLIterativeLeastSquare3DNewton.c

OBJS := ${SRCS:%.c=$(libdir)/%.o}
INCS := ${SRCS:%.c=$(incdir)/%.h}
PKGINCS := ${SRCS:%.c=$(pkgincdir)/kagali/%.h}

CFLAGS  := -ggdb -std=c99 -Wall -O2 -I$(pkgincdir)
LDFLSGS := -ggdb -std=c99 -Wall -O2

MKPROTO  := $(scriptdir)/mkproto.pl
KGLVCSID := $(scriptdir)/kglvcsid.pl
GITDIR   := $(pkgdir)/.git

.PRECIOUS: $(INCS)

all: $(OBJS)

$(incdir)/%.h : %.c %.h.in
	$(MKPROTO) $< $(incdir)
	cd $(incdir) && ln -sf `pwd`/`basename $@` $(pkgincdir)/kagali

$(libdir)/%.o : %.c $(incdir)/%.h
	gcc $(CFLAGS) -I/opt/local/include -c -o$@ $<

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(INCS) *~
	rm -f $(PKGINCS)
