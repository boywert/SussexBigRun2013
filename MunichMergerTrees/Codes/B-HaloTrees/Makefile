
EXEC   = B-HaloTrees

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o

INCL   = allvars.h  proto.h  Makefile

#OPTIONS  += -DLONGIDS
#OPTIONS  += -DSAVE_MASS_TAB         # writes MassTab of individual species into output

SYSTYPE="tiberius"
#SYSTYPE="cosma-iv"
#SYSTYPE="Octor"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"


CC       =    cc        # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall  # optimization and warning flags (default)


ifeq ($(SYSTYPE),"tiberius")
CC = cc
OPTIMIZE     = -O3
endif


ifeq ($(SYSTYPE),"cosma-iv")
CC = icc
OPTIMIZE     = -O3 -Wall -shared-intel -mcmodel=medium
endif

ifeq ($(SYSTYPE),"Octor")
CC = gcc
OPTIMIZE = -m64 -O3 -Wall
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

### uncomment PARSING for using support of ProC parameter files ###
#PARSING="PROCCOMPATIBLE"

ifeq ($(PARSING), "PROCCOMPATIBLE")
OBJS += paramparser.o
OPTIONS += -D PROCCOMPATIBLE
endif

###

LIBS   =   -lm -g

CFLAGS =   -g $(OPTIONS)  $(OPTIMIZE) 

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS)  $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)







