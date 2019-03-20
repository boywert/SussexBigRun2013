EXEC   = B-BaseTree

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o

INCL   = allvars.h  proto.h  Makefile


OPTIONS += -DFLAG_GROUP_VELDISP
#OPTIONS += -DLONGIDS
#OPTIONS += -DIDS_HAVE_GAPS            
OPTIONS += -DHAVE_HDF5 
#OPTIONS += -DH5_USE_16_API
#OPTIONS += -DHAVE_NUMPART_HIGHWORD
#OPTIONS += -DSKIP_BY_WEIGHT           # Do A->C assigning also in case that substructure isassigned to be merged to larger structure in B
#OPTIONS += -DALPHA=0.8                # Exponent in weighting function
#OPTIONS += -DBACKWARD_CHECKING        # Check also backwards (B->A) dependencies and does reassignement of decendent in case of halos appearing 
#OPTIONS += -DHALO_SIZE_INCREASE_FOR_SWITCHING=1.5    # Howmuch bigger a "newborn" Halo has to be to redirect in case of detected confusion

SYSTYPE="tiberius"
#SYSTYPE="pollux"
#SYSTYPE="cosma-iv"
#SYSTYPE="Octor"
#SYSTYPE="VIP"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"


CC       =    cc        # sets the C-compiler (default)
OPTIMIZE =   -g -Wall    # optimization and warning flags (default)
HDF5_LIBS    =  -lhdf5
HDF5_LDFLAGS = -L/pactusers/cs390/local/lib
HDF5_INCLUDE = -I/pactusers/cs390/local/include

ifeq ($(SYSTYPE),"tiberius")
CC = gcc
OPTIMIZE     = -O3 -Wall 
HDF5_LIBS    = -lhdf5
HDF5_INCLUDE = -I/users/home/boyd/local/include -DH5_USE_16_API
HDF5_LDFLAGS  =  -L/users/home/boyd/local/lib  -lhdf5
#HDF5_INCLUDE = #-I/pactusers/cs390/local/include
endif



ifeq ($(SYSTYPE),"pollux")
CC = icc
OPTIMIZE     = -O3 -Wall -shared-intel -mcmodel=medium -openmp     
#OPTIMIZE     = -O0 -g -Wall -shared-intel -mcmodel=medium #-openmp                                               
HDF5_LIBS    = -lhdf5
HDF5_LDFLAGS = $(LDFLAGS)
HDF5_INCLUDE = $(CPPFLAGS)
endif

ifeq ($(SYSTYPE),"cosma-iv")
CC = icc
OPTIMIZE     = -O3 -Wall #-openmp
HDF5_LIBS    = -lhdf5
HDF5_LDFLAGS = $(LDFLAGS)
HDF5_INCLUDE = $(CPPFLAGS)
endif

ifeq ($(SYSTYPE),"Octor")
CC = gcc
OPTIMIZE = -m64 -O3 -Wall
endif

ifeq ($(SYSTYPE),"VIP")
CC       =   cc_r
OPTIMIZE =   -q64  -qsmp #  -O5 -qstrict -qipa 
endif


ifeq ($(SYSTYPE),"Regatta")
CC       =   cc_r
OPTIMIZE =   -q64 -qsmp #  -O5 -qstrict -qipa 
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


LIBS   =   -lm -g $(HDF5_LIBS) $(HDF5_LDFLAGS)

CFLAGS =   -g $(OPTIONS)  $(OPTIMIZE) $(HDF5_INCLUDE)
CFLAGS += $(LIBS)
$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS)  $(CFLAGS)  -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)









