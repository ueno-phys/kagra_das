CC := gcc
pkgdir := ../../..
pkgincdir := $(pkgdir)/include
pkglibdir := $(pkgdir)/lib
scriptdir := $(pkgdir)/scripts
incdir := ../include

libdir := ../lib

SRCS := \
	KGLError.c \
	KGLVersion.c \
	nrutil.c

OBJS := ${SRCS:%.c=$(libdir)/%.o}

CFLAGS  := -ggdb -std=c99 -Wall -O2 -I$(pkgincdir)
LDFLSGS := -ggdb -std=c99 -Wall -O2

MKPROTO  := $(scriptdir)/mkproto.pl
KGLVCSID := $(scriptdir)/kglvcsid.pl
GITDIR   := $(pkgdir)/.git

all: $(OBJS)

$(incdir)/%.h : %.c %.h.in
	$(MKPROTO) $< $(incdir)
	cd $(incdir) && ln -sf `pwd`/`basename $@` $(pkgincdir)/kagali

KGLVersion.h: $(GITDIR)/index
	$(KGLVCSID) KGL

$(libdir)/KGLError.o: KGLError.c $(incdir)/KGLError.h KGLVersion.h configStd.h
	gcc $(CFLAGS) -c -o$@ $<

$(libdir)/KGLVersion.o: KGLVersion.c KGLVersion.h
	gcc $(CFLAGS) -c -o$@ $<

$(libdir)/nrutil.o: nrutil.c $(incdir)/nrutil.h
	gcc $(CFLAGS) -c -o$@ $<

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(incdir)/KGLError.h KGLVersion.h *~
	rm -f $(pkgincdir)/kagali/KGLError.h
