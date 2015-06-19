# This is the Makefile for the library routines used in for the random
# number generator

# setenv MallocHelp   -- to see help about tracking memory

OS ?= macosx
# OS ?= linux

CC := gcc
CCPP := g++

INCLUDE := -I.
LIBS := -L.

LIBM := -lm
LIBGSL := -lgsl -lgslcblas
ifeq ($(OS),macosx)
#  the first one protects the -framework vecLib: only needed for xlf/xlc
  OSDEFS := -DMACOSX
  LIBS := -L/sw/lib $(LIBS)
  INCLUDE := -I/sw/include $(INCLUDE)
endif

CFLAGS ?= -O5
# have_inline defined for gsl to use extern inline functions when possible:
CXXFLAGS ?= -O5 -DHAVE_INLINE $(OSDEFS)

TARGET = wsflux

INCLUDES = *.H
# INCLUDES = io.H multipoly.H eam.H

.PHONY: all clean archive

all: $(TARGET)

wslux:	wsflux.C $(INCLUDES)

.SUFFIXES: .C .c .o

.C: $< $(INCLUDES)
	$(CCPP) -Wconversion $(CXXFLAGS) $(INCLUDE) $< -o $@ $(LIBS)
#	$(CCPP) -Wconversion $(CXXFLAGS) $(INCLUDE) $< -o $@ $(LIBS) $(LIBGSL) $(LIBM)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.C.o:
	$(CCPP) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm -f *.o *.a

SOURCES = *.C *.H
MAKEFILES = Makefile
ARCHIVEFILES = $(SOURCES) $(MAKEFILES)

archive:
	@makearchive -z $(ARCHIVEFILES)
