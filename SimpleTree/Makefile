CC =  mpicc
CFLAGS = -g -Wall -openmp

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
#OPT += -DTOPLEVELONLY #use only toplevel halos
VARS = 
VARS += -Dnpart_box 5159780352
VARS += -Dboxsize 47000.0
VARS += -Ddomain_per_dim 6
VARS += -DINPUTDIR "/mnt/lustre/scratch/cs390/Sussexbigrun2013/cubepm_130315_6_1728_47Mpc_ext2/results"
VARS += -DOUTPUTDIR "/mnt/lustre/scratch/cs390/Sussexbigrun2013/cubepm_130315_6_1728_47Mpc_ext2/mergertree"
CFLAGS += $(OPT)



all: maketree

maketree: $(SRCS) $(HEADERS) $(LIBSRC)
	$(CC) $(CFLAGS) $(SRCS) $(HEADERS) $(LIBSRC) -o maketree

clean:
	rm -f maketree
