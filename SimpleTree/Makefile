CC =  gcc  -I/opt/mpi/bullxmpi/1.1.16.5/include -pthread -L/opt/mpi/bullxmpi/1.1.16.5/lib -lmpi -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
CFLAGS = -g -lm -Wall -fopenmp

SRCDIR = src

HEADERS = 

LIBSRC = $(SRCDIR)/libio/common_io.h $(SRCDIR)/libio/common_io.c
LIBSRC += $(SRCDIR)/libio/readsussexbigrun.h $(SRCDIR)/libio/readsussexbigrun.c
LIBSRC += $(SRCDIR)/libmemmgr/memmgr.h $(SRCDIR)/libmemmgr/memmgr.c
LIBSRC += $(SRCDIR)/libhash/hash.h $(SRCDIR)/libhash/hash.c

SRCS = $(SRCDIR)/main.c $(SRCDIR)/main.h
SRCS += $(SRCDIR)/common.c $(SRCDIR)/common.h
SRCS += $(SRCDIR)/treeformat.c $(SRCDIR)/treeformat.h
SRCS += $(SRCDIR)/cosmology.h $(SRCDIR)/cosmology.c
SRCS += $(SRCDIR)/descendant.h $(SRCDIR)/descendant.c
SRCS += $(SRCDIR)/MPI_distribute.h $(SRCDIR)/MPI_distribute.c

OPT = 
OPT += -DCUBEP3M #enable domains in cubep3m

CFLAGS += $(OPT)

all: maketree

maketree: $(SRCS) $(HEADERS) $(LIBSRC)
	$(CC) $(CFLAGS) $(SRCS) $(HEADERS) $(LIBSRC) -o maketree

clean:
	rm -f maketree
