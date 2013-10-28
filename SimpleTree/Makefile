CC =  mpicc
CFLAGS = -g -openmp

SRCDIR = src
BINDIR = bin
HEADERS = 

LIBSRC = $(SRCDIR)/libio/common_io.h $(SRCDIR)/libio/common_io.c
LIBSRC += $(SRCDIR)/libio/readsussexbigrun.h $(SRCDIR)/libio/readsussexbigrun.c
LIBSRC += $(SRCDIR)/libmemmgr/memmgr.h $(SRCDIR)/libmemmgr/memmgr.c
LIBSRC += $(SRCDIR)/libhash/hash.h $(SRCDIR)/libhash/hash.c

MAINSRCS = $(SRCDIR)/main.c $(SRCDIR)/main.h 
TESTSRCS = $(SRCDIR)/test.c $(SRCDIR)/main.h 
CCCSRCS = $(SRCDIR)/chunk2cubep3m.c $(SRCDIR)/main.h 

SRCS = 
SRCS += $(SRCDIR)/readconfig.c $(SRCDIR)/readconfig.h
SRCS += $(SRCDIR)/common.c $(SRCDIR)/common.h
SRCS += $(SRCDIR)/treeformat.c $(SRCDIR)/treeformat.h
SRCS += $(SRCDIR)/cosmology.h $(SRCDIR)/cosmology.c
SRCS += $(SRCDIR)/descendant.h $(SRCDIR)/descendant.c
SRCS += $(SRCDIR)/MPI_distribute.h $(SRCDIR)/MPI_distribute.c

MAINSRCS += $(SRCS)
TESTSRCS += $(SRCS)
CCCSRCS += $(SRCS)
OPT = 
OPT += -DCUBEP3M #enable domains in cubep3m
#OPT += -DTOPLEVELONLY #use only toplevel halos

CFLAGS += $(OPT)



all: maketree test chunk2cubep3m

maketree: $(MAINSRCS) $(HEADERS) $(LIBSRC)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(MAINSRCS) $(HEADERS) $(LIBSRC) -o $(BINDIR)/maketree

test: $(TESTSRCS) $(HEADERS) $(LIBSRC)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(TESTSRCS) $(HEADERS) $(LIBSRC) -o $(BINDIR)/test

chunk2cubep3m: $(CCCSRCS) $(HEADERS) $(LIBSRC)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(TESTSRCS) $(HEADERS) $(LIBSRC) -o $(BINDIR)/chunk2cubep3m
clean:
	rm -f $(BINDIR)/*
