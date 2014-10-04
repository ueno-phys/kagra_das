staticlib := lib/libkagali.a

subdirs := \
	kagali/std \
	kagali/fft \
	kagali/detector \
	kagali/cbc \
	kagali/nha \
	kagali/clock

objdirs := $(subdirs:%=%/lib)
objs := $(shell find $(objdirs) -name '*.o')

all: all-recursive
	make $(staticlib)

all-recursive:
	for dir in $(subdirs); do make -C $$dir/src || exit 1 ; done

$(staticlib): $(objs)
	ar cru $@ $?
	ranlib $@

clean: clean-recursive
	rm -f $(staticlib)

clean-recursive:
	for dir in $(subdirs); do make -C $$dir/src clean; done

distclean: distclean-recursive
	rm -f $(staticlib) *~

distclean-recursive:
	for dir in $(subdirs); do make -C $$dir/src distclean; done
