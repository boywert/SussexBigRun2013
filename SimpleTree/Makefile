CC = gcc
CFLAGS = -lm -fopenmp
SRCDIR = src
HEADERS = $(SRCDIR)/main.h $(SRCDIR)/treeformat.h
LIBSRC = $(SRCDIR)/libio/common_io.h $(SRCDIR)/libio/common_io.c
LIBSRC += $(SRCDIR)/libio/readsussexbigrun.h $(SRCDIR)/libio/readsussexbigrun.c
LIBSRC += $(SRCDIR)/libmemmgr/memmgr.h $(SRCDIR)/libmemmgr/memmgr.c
SRCS = $(SRCDIR)/main.c 

all: maketree

maketree: $(SRCS) $(HEADERS) $(LIBSRC)
	$(CC) $(CFLAGS) $(SRCS) $(HEADERS) $(LIBSRC) -o maketree