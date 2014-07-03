EXEC   = L-Galaxies

OBJS   = ./code/main.o ./code/io_tree.o ./code/read_parameters.o ./code/allvars.o 
	 


INCL   = ./code/allvars.h  ./code/proto.h  Makefile


#OPT     += -DMR
# OPT     += -DMRII # MR + MRII

OPT += -DTREES
OPT += -DTREE_IDS
OPT += -DAUX


#OPT += -DHALO_MODEL

# SYSTYPE = "MPA"
# SYSTYPE = "MPA-AMD"
# SYSTYPE = "Regatta"
# SYSTYPE = "RZG_LinuxCluster"
# SYSTYPE = "RZG_LinuxCluster-gcc"
# SYSTYPE = "Solaris"
# SYSTYPE = "OPA-Cluster64"
# SYSTYPE = "cygwin"
# SYSTYPE = "COSM"

CC       = gcc				# sets the C-compiler (default)
CC_MD       = gcc				# sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall 	# optimization and warning flags (default)

ifeq ($(SYSTYPE),"COSM")
#add this in your .cshrc
#module load gsl/1.14
#module load fftw/2.1.5
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =  mpicc
else
CC       =  icc
endif
OPTIMIZE =  -O3 -Wall -m64
GSL_INCL =  -I$(GSL_HOME)/include
GSL_LIBS =  -L$(GSL_HOME)/lib -Xlinker -R -Xlinker $(GSL_HOME)/lib 
FFTW_INCL=  -I$(FFTW_HOME)/include
FFTW_LIBS=  -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =  
endif


ifeq ($(SYSTYPE),"MPA-AMD")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       = /opt/openmpi-gnu-1.2.6/bin/mpicc  -m64  -L/opt/openmpi-gnu-1.2.6/lib  
else
CC       =   cc -pg
#CC       =   cc 
endif
endif

ifeq ($(SYSTYPE),"OPA-Cluster64")
CC       =  mpiccg
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCCg -m64
endif
OPT     +=  -DNOCALLSOFSYSTEM
OPTIMIZE =  -O3 -Wall -m64
GSL_INCL =  -I/afs/rzg/bc-b/vrs/opteron64/include
GSL_LIBS =  -L/afs/rzg/bc-b/vrs/opteron64/lib  -Wl,"-R /afs/rzg/bc-b/vrs/opteron64/lib"
FFTW_INCL=  -I/afs/rzg/bc-b/vrs/opteron64/include
FFTW_LIBS=  -L/afs/rzg/bc-b/vrs/opteron64/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =   mpicc
else
CC       =   cc
endif
GSL_INCL = -I/usr/common/pdsoft/include
GSL_LIBS = -L/usr/common/pdsoft/lib -static
endif


ifeq ($(SYSTYPE),"Opteron")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =   mpicc
else
CC       =   cc
endif
GSL_INCL = -I/usr/local/include
GSL_LIBS = -L/usr/local/lib -static
endif



ifeq ($(SYSTYPE),"Regatta")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =   mpcc_r
else
CC       =   cc_r
endif
OPTIMIZE =   -g -q64 -qflttrap=zerodivide:invalid:enable  #  -O5 -qstrict -qipa 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include  -I/afs/rzg/@sys/lib/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib                
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =   mpicci
else
CC       =   icc
endif
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/lib                -static
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
ifeq (PARALLEL,$(findstring PARALLEL,$(OPT)))
CC       =   mpiccg
else
CC       =   gcc
endif
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/lib
endif

ifeq ($(SYSTYPE),"Solaris")
CC       =   cc   # sets the C-compiler
OPTIMIZE =   -i -fast -xvector -xarch=v9b -xchip=ultra3 -xcache=64/32/4:8192/512/1 -I/opt/local/include
endif

ifeq ($(SYSTYPE),"cygwin")
CC       =   gcc   # sets the C-compiler
GSL_INCL = -I/usr/local/include
GSL_LIBS = -L/usr/local/lib -static
#GSL_INCL = -I/usr/include
#GSL_LIBS = -L/usr/lib -static
endif


LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =   -g $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E ./code/allvars.h -o ./code/allvars.i
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/GALAXY_OUTPUT_2_TypeString.awk > ./awk/L-Galaxies_Types.txt
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/GALAXY_OUTPUT_2_DDL.awk > ./awk/L-Galaxies_DDL.sql
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_struct.awk > /afs/mpa/data/bmh20/SAM/LGalaxy.pro
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_hists.awk >> /afs/mpa/data/bmh20/SAM/LGalaxy.pro


	
