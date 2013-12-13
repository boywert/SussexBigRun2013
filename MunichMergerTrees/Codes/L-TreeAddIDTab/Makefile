EXEC   = L-TreeAddIDTab

OBJS   = main.o io_tree.o allvars.o mymalloc.o read_parameters.o

INCL   = allvars.h proto.h  Makefile

#OPTIONS   +=  -DSAVE_MASS_TAB

#SYSTYPE="cosma-iv"
#SYSTYPE="Octor"
#SYSTYPE="Rosalind"
#SYSTYPE="Archi"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"
SYSTYPE="tiberius"

CC       =    cc -g       # sets the C-compiler (default)
OPTIMIZE =    -Wall    # optimization and warning flags (default)

ifeq ($(SYSTYPE),"cosma-iv")
CC = icc
OPTIMIZE     = -O3 -Wall -shared-intel -mcmodel=medium
endif

ifeq ($(SYSTYPE),"tiberius")
CC = cc
OPTIMIZE     = -O3 -Wall 
endif


ifeq ($(SYSTYPE),"Octor")
CC = gcc
OPTIMIZE = -m64 -O3 -Wall
endif

ifeq ($(SYSTYPE),"Rosalind")
CC = cc
OPTIMIZE = -m64 -O3
#OPTIMIZE = -fast -xdepend -m64 -xarch=sse2a
endif

ifeq ($(SYSTYPE),"Archi")
CC = cc
OPTIMIZE =  -O3 -Wall
endif

ifeq ($(SYSTYPE),"Regatta")
CC       =   cc
OPTIMIZE =   -g -q64 #  -O5 -qstrict -qipa 
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   icc
OPTIMIZE =   -O3 
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   gcc
OPTIMIZE =   -O3 
endif

ifeq ($(SYSTYPE),"Solaris")
CC       =   cc   # sets the C-compiler
OPTIMIZE =   -i -fast -xvector -xarch=v9b -xchip=ultra3 -xcache=64/32/4:8192/512/1 -I/opt/local/include
endif


LIBS   =   -lm -g

CFLAGS =   -g $(OPTIONS)  $(OPTIMIZE)  

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
