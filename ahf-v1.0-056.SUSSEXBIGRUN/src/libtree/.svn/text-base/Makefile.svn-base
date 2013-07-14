# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS+=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# If a different MAKE is specified use that one, otherwise try your luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

OBJS = generate_tree.o

# Catch empty make call
all:
	${MAKE} libtree.a

# Cleaning up
clean:
	rm -f libtree.a $(OBJS)

# The library itself
libtree.a: $(OBJS)
	${AR} -r libtree.a $(OBJS)

# The individual routines stuff
generate_tree.o:	generate_tree.c
	$(CC) $(CFLAGS) -c -o generate_tree.o generate_tree.c

