EXEC   = L-TreeAddPosTab

OBJS   = main.o  allvars.o mymalloc.o read_parameters.o hdf5_io.o

INCL   = allvars.h proto.h  Makefile


#OPTIONS += -DLONGIDS
#OPTIONS += -DHAVE_HDF5 
#OPTIONS += -DH5_USE_16_API
#OPTIONS += -DPATCH_IDS


#SYSTYPE="pollux"
SYSTYPE="tiberius"
#SYSTYPE="Octor"
#SYSTYPE="Rosalind"
#SYSTYPE="Archi"
#SYSTYPE="Opteron"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"


CC       =    cc -g       # sets the C-compiler (default)
OPTIMIZE =    -Wall    # optimization and warning flags (default)
HDF5_LIBS =
HDF5_LDFLAGS =
HDF5_INCLUDE =

ifeq ($(SYSTYPE),"tiberius")
CC = gcc
OPTIMIZE     = -O3 -Wall 
HDF5_LIBS    = #-lhdf5
HDF5_LDFLAGS = #-L/pactusers/cs390/local/lib
HDF5_INCLUDE = #-I/pactusers/cs390/local/include
endif



ifeq ($(SYSTYPE),"pollux")
CC = icc
OPTIMIZE     = -O3 -Wall -shared-intel -mcmodel=medium  
HDF5_LIBS    = -lhdf5
HDF5_LDFLAGS = $(LDFLAGS)
HDF5_INCLUDE = $(CPPFLAGS)
endif

ifeq ($(SYSTYPE),"cosma-iv")
CC = icc
OPTIMIZE     = -O3 -Wall
HDF5_LIBS    = -lhdf5
HDF5_LDFLAGS = $(LDFLAGS)
HDF5_INCLUDE = $(CPPFLAGS)
endif

ifeq ($(SYSTYPE),"Octor")
CC = cc
OPTIMIZE = -m64 -O3
endif

ifeq ($(SYSTYPE),"Rosalind")
CC = cc
OPTIMIZE = -m64 -O3
#OPTIMIZE = -fast -xdepend -m64 -xarch=sse2a
endif

ifeq ($(SYSTYPE),"Archi")
CC = cc
OPTIMIZE =  -O3
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

ifeq ($(SYSTYPE),"Opteron")
CC       =   gcc
OPTIMIZE =   -W -Wall -O2 -m64 
endif

ifeq ($(SYSTYPE),"Solaris")
CC       =   cc   # sets the C-compiler
OPTIMIZE =   -i -fast -xvector -xarch=v9b -xchip=ultra3 -xcache=64/32/4:8192/512/1 -I/opt/local/include
endif


LIBS   =   -g $(HDF5_LIBS) -lm

CFLAGS =    $(OPTIONS)  $(OPTIMIZE) $(HDF5_LDFLAGS) $(HDF5_INCLUDE)  

CFLAGS += $(LIBS)
$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS)  -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
