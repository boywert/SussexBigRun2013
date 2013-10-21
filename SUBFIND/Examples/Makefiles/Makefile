#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=Gadget_new
# 
#-----------------------------------------------------------------

#--------------------------------------- Select target Computer

#SYSTYPE="Darwin"
#SYSTYPE="odyssey"
#SYSTYPE="OpenSuse"
#SYSTYPE="OpenSuse64"
#SYSTYPE="Genius"
#SYSTYPE="MBM"
#SYSTYPE="Magny"
#SYSTYPE="Solaris_amd64"
#SYSTYPE="HLRB2"
#SYSTYPE="SUNNYVALE"
#SYSTYPE="McKenzie"
#SYSTYPE="CITA"
#SYSTYPE="Stella"
#SYSTYPE="Sauron"
#SYSTYPE="Sauron-gcc"
#SYSTYPE="Mako"
#SYSTYPE="MPA"
#SYSTYPE="Regatta"
#SYSTYPE="VIP"
#SYSTYPE="SARA"
#SYSTYPE="CINECA_SP6"
#SYSTYPE="Ingeld_LinuxCluster"
#SYSTYPE="Zijl"
#SYSTYPE="hpcf"
#SYSTYPE="opus"
#SYSTYPE="OpteronMPA-Gnu"
#SYSTYPE="OpteronMPA-Intel"
#SYSTYPE="OPA-Cluster64-Intel"
#SYSTYPE="Warp"
#SYSTYPE="Kraken"
#SYSTYPE="Kraken_FFTW3"
#SYSTYPE="Octopus"
#SYSTYPE="Ranger_pgi"
#SYSTYPE="Ranger_intel"
#SYSTYPE="cosma-iv"
SYSTYPE="tiberius"
#----------------------------------------------------------------------------------------------

CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat gadgetconfig.h)


CC       = mpicc        # sets the C-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =   mpiCC     # default C++ compiler
OPTIMIZE =   -g 
OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems
endif

ifeq (VORONOI,$(findstring VORONOI,$(CONFIGVARS)))
GMPLIB   =  -lgmp
endif

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES =  -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES =  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES =  -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif



#----------------------------------------------------------------------------------------------

ifeq ($(SYSTYPE),"cosma-iv")
CC	 =  /pactusers/cs390/mpirun/mpicc
OPTIMIZE =  #-O3 -xHOST #-g -openmp -Wall
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
MPICHLIB = 
endif


ifeq ($(SYSTYPE),"cosma-iv")
CC       =  mpicc
OPTIMIZE =  -O3 -xHOST #-g -openmp -Wall
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
MPICHLIB = 
endif

ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc -g -O2 -xW -ipo -Wall
OPTIMIZE =
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc -g -O2 -tp barcelona-64
OPTIMIZE =
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
endif



ifeq ($(SYSTYPE),"Octopus")
CC       =  mpicc
OPTIMIZE = -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS= -lfftw3
MPICHLIB = -lmpi
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc 
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC  
endif
#gnu
OPT      +=  -DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
OPTIMIZE =   -O3 -g -Wall -m64
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -I/n/sw/hdf5-1.8.5_gnu/include
HDF5LIB  =  -L/n/sw/hdf5-1.8.5_gnu/lib -lhdf5
#intel
#OPT      +=  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
#OPTIMIZE =   -parallel -O3 -ipo -funroll-loops -Wall
#GSL_INCL =
#GSL_LIBS =
#FFTW_INCL=
#FFTW_LIBS=
#MPICHLIB =
#HDF5INCL =
#HDF5LIB  =  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Darwin")
CC       =  mpicc -g -O2 -Wall
OPTIMIZE =  
GSL_INCL = -I/sw/include 
GSL_LIBS = -L/sw/lib 
FFTW_INCL= -I/sw/include 
FFTW_LIBS= -L/sw/lib 
MPICHLIB = -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"Genius")
CC       =  mpixlc_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/u/vrs/BlueGene/Libs/include
GSL_LIBS =  -L/u/vrs/BlueGene/Libs/lib -static
FFTW_INCL=  -I/u/vrs/BlueGene/Libs/include
FFTW_LIBS=  -L/u/vrs/BlueGene/Libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Solaris_amd64")
CC       = /opt/SUNWhpc/bin/mpicc
OPTIMIZE = -xO4 -xchip=native -xarch=amd64a
#OPTIMIZE =   -g -xarch=amd64a
GSL_INCL = $(GSL64INCL)
GSL_LIBS = $(GSL64LIB)
FFTW_INCL= $(FFTW64INCL)
FFTW_LIBS= $(FFTW64LIB)
FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
R_PATH   = $(R64LIB)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)  $(GSL_LIBS)
HDF5INCL =
HDF5LIB  =
LIBS     =   -lm  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) $(R64LIB)
endif


ifeq ($(SYSTYPE),"HLRB2")
CC       =  icc   # sets the C-compiler
OPTIMIZE =  -g -O3 -ipo  -mtune=itanium2 -mcpu=itanium2 -mp
GSL_INCL = -I/home/hlrb2/h1091/h1091aa/Libs/include
GSL_LIBS = -L/home/hlrb2/h1091/h1091aa/Libs/lib -Xlinker -R -Xlinker /home/hlrb2/h1091/h1091aa/Libs/lib
FFTW_INCL= -I/home/hlrb2/h1091/h1091aa/Libs/include
FFTW_LIBS= -L/home/hlrb2/h1091/h1091aa/Libs/lib -Xlinker -R -Xlinker /home/hlrb2/h1091/h1091aa/Libs/lib
MPICHLIB = -lmpi
HDF5INCL = 
HDF5LIB  = 
endif


ifeq ($(SYSTYPE),"MBM") 
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/sw/tap/include
GSL_LIBS =  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
FFTW_INCL=  -I/sw/tap/include
FFTW_LIBS=  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif

# modules for Magny
# module add mvapich2/gcc/64/1.4.1-qlc

ifeq ($(SYSTYPE),"Magny") 
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/hits/tap/sw/libs/include
GSL_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
FFTW_INCL=  -I/hits/tap/sw/libs/include
FFTW_LIBS=  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
GMP_INCL =  -I/hits/tap/sw/libs/include
GMP_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
MPICHLIB =
HDF5INCL =  -I/hits/tap/sw/libs/include -DH5_USE_16_API
HDF5LIB  =  -L/hits/tap/sw/libs/lib -lhdf5 
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"OpenSuse")
CC       =  mpicc
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"OpenSuse64")
CC       =  mpicc
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib64/mpi/gcc/openmpi/lib64 -Xlinker -R -Xlinker /usr/lib64/mpi/gcc/openmpi/lib64 -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"SUNNYVALE")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib64
FFTW_INCL= -I${MCKENZIE_FFTW_INC_PATH}
FFTW_LIBS= -L${MCKENZIE_FFTW_LIB_PATH}
MPICHLIB =  -L/opt/lam-7.1.2-intel/lib -lmpi
HDF5INCL = -I${MCKENZIE_HDF_INC_PATH}
HDF5LIB  = -L${MCKENZIE_HDF_LIB_PATH} -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"McKenzie")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib
FFTW_INCL= -I/opt/fftw/intel_8.1/2.1.5/include
FFTW_LIBS= -L/opt/fftw/intel_8.1/2.1.5/lib #-ldrfftw_mpi
MPICHLIB = -L/opt/lam-7.1.2b24-g77/lib -lmpi
HDF5INCL = -I/opt/hdf5-oscar-1.6.4/include
HDF5LIB  = -L/opt/hdf5-oscar-1.6.4/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif


# modules for OPA  
# module load intel
# module load impi
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5-serial

ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  
CC       =   mpiicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64 -wd981 -wd2259 -wd1572 -ipo
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR) -Xlinker -R -Xlinker $(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib -Xlinker -R -Xlinker $(FFTW_HOME)/lib
MPICHLIB =  
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -Xlinker -R -Xlinker $(HDF5_HOME)/lib -lhdf5 -lz
#OPT      +=  -DOLD_HDF5
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"opus")
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/snfs1/cgwork/dgad0043/include
GSL_LIBS =  -L/snfs1/cgwork/dgad0043/lib
FFTW_INCL=  -I/snfs1/cgwork/dgad0043/include
FFTW_LIBS=  -L/snfs1/cgwork/dgad0043/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC  -m64
endif
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"OpteronMPA-Intel")
CC       =  /usr/local/appl/openmpi-intel/bin/mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  /usr/local/appl/openmpi-intel/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"OpteronMPA-Gnu")
CC       =  /usr/local/appl/openmpi-gnu/bin/mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  /usr/local/appl/openmpi-gnu/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -g -Wall -fopenmp 
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"Mako")
CC       =  mpicc   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -O3 -march=athlon-mp  -mfpmath=sse
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Regatta")
CC       =  mpcc_r -g # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPTIMIZE =  -qstrict -q64 -qcpluscmt  #-O3 -qipa
GSL_INCL = -I/afs/rzg/u/vrs/gsl_psi64/include
GSL_LIBS = -L/afs/rzg/u/vrs/gsl_psi64/lib
FFTW_INCL= -I/afs/rzg/u/vrs/fftw_psi64/include
FFTW_LIBS= -L/afs/rzg/u/vrs/fftw_psi64/lib  -q64 # -qipa
MPICHLIB =
HDF5INCL = -I/afs/rzg/u/vrs/hdf5_psi64/include
HDF5LIB  = -L/afs/rzg/u/vrs/hdf5_psi64/lib  -lhdf5 -lz
endif

# modules for VIP
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5

ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g -qcpluscmt # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline -qstrict -qcpluscmt -qipa # -qhot -qipa=level=2 -qfloat=hdflt
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
OPT     +=   -DPOWER6 -DALT_QSORT -DPAD_STRUCTURES
endif

# modules for SARA
# module load c/ibm/11.1
# module load fftw2/2.1.5
# module load gsl/1.11
# module load hdf5

ifeq ($(SYSTYPE),"SARA")
CC       =  mpcc -g -qcpluscmt -qsmp -qreport # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline -qstrict -qhot -qipa=level=2 -qfloat=hsflt
GSL_INCL = -I$(SARA_GSL_INCLUDE)
GSL_LIBS = -L$(SARA_GSL_LIB)
FFTW_INCL= -I$(SARA_FFTW2_INCLUDE)
FFTW_LIBS= -L$(SARA_FFTW2_LIB)
MPICHLIB =
HDF5INCL = -I$(SARA_HDF5_INCLUDE)
HDF5LIB  = -L$(SARA_HDF5_LIB) -lhdf5 -lz -static
OPT     +=   -DPOWER6 -DALT_QSORT -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"CINECA_SP6")
CC       =   mpcc # sets the gcc C-compiler
#to be used for deugging: -g -qcheck=all -qdbxextra -qheapdebug -qkeepparm -qoptdebug
ifndef DEBUG
OPTIMIZE =  -O3 -q64  -qhot -qstrict -qfloat=hsflt -qarch=pwr6 -qtune=pwr6 -qcpluscmt  #-qipa=level=2:inline=noauto:list
else
OPTIMIZE = $(DEBUG) -qarch=pwr6 -qtune=pwr6 -qcpluscmt
endif

GSL_INCL = -I/cineca/prod/libraries/gsl/1.9/xl--10.1/include
GSL_LIBS = -L/cineca/prod/libraries/gsl/1.9/xl--10.1/lib
FFTW_INCL= -I/cineca/prod/libraries/fftw/2.1.5/xl--10.1/include
FFTW_LIBS= -L/cineca/prod/libraries/fftw/2.1.5/xl--10.1/lib
HDF5INCL = -I/sp6/userexternal/sborgani/code/HDF5/include
HDF5LIB  = -L/sp6/userexternal/sborgani/code/HDF5/lib -lhdf5
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Ingeld_LinuxCluster")
CC       =   mpicc  # sets the C-compiler
OPTIMIZE =   -O3 -Wall

GSL_INCL = -I/home/patricia/my-gsl/include
GSL_LIBS = -L/home/patricia/my-gsl/lib -static
FFTW_INCL= -I/home/patricia/my-fftw/include
FFTW_LIBS= -L/home/patricia/my-fftw/lib
endif


ifeq ($(SYSTYPE),"hpcf")
CC       =   mpicc  # sets the C-compiler
OPT     +=  -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPTIMIZE =  -O3
GSL_INCL = -I/home/gadget/Libs/include
GSL_LIBS = -L/home/gadget/Libs/lib
FFTW_INCL= -I/home/gadget/Libs/include
FFTW_LIBS= -L/home/gadget/Libs/lib
endif


ifeq ($(SYSTYPE),"Warp")
CC       =   /opt/intel/impi/3.1/bin64/mpiicc
OPTIMIZE =  -xW
GSL_INCL = -I/opt/gsl/impi/include/gsl
GSL_LIBS = -L/opt/gsl/impi/lib64 -Wl,-rpath,/opt/gsl/impi/lib64
FFTW_INCL= -I/opt/fftw/impi_intel-3.1/include
FFTW_LIBS= -L/opt/fftw/impi_intel-3.1/lib64 -Wl,-rpath,/opt/fftw/impi_intel-3.1/lib64
MPICHLIB = -L/opt/intel/impi/3.1/lib64 -Wl,-rpath,/opt/intel/impi/3.1/lib64
endif


ifeq ($(SYSTYPE),"Kraken")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/opt/fftw/2.1.5/cnos/include
FFTW_LIBS=  -L/opt/fftw/2.1.5/cnos/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
endif

ifeq ($(SYSTYPE),"Kraken_FFTW3")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/include
FFTW_LIBS=  -L/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))    # fftw installed with type prefix?
  FFTW_LIBNAMES = $(FFTW_LIBS) -lfftw3_mpi -lfftw3_threads -lfftw3
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))
  FFTW_LIBNAMES =  -lfftw3_mpi -lfftw3_threads -lfftw3
else
  FFTW_LIBNAMES =  -lfftw3f_mpi -lfftw3f_threads -lfftw3f
endif
endif

endif



ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif


ifeq (CHEMCOOL,$(findstring CHEMCOOL,$(CONFIGVARS)))
FOPT	+=  -WF,-DCHEMCOOL
OPT	+=  -DCHEMISTRYNETWORK=1
FOPT	+=  -WF,-DCHEMISTRYNETWORK=1
endif


ifneq (MODGRAV,$(findstring MODGRAV,$(CONFIGVARS))) # standard gravity object files
GRAVITY_OBJS  =  forcetree.o darkenergy.o pm_periodic.o
else # replacements used in the modified gravity model 
GRAVITY_OBJS  =  modgrav/forcetree.o modgrav/darkenergy.o modgrav/pm_periodic.o modgrav/modgrav_background.o
endif

L3_OBJS = 

ifeq (SORT_FROM_L3,$(findstring SORT_FROM_L3,$(CONFIGVARS))) # add L3 parallel sorting
L3_OBJS  +=  ImportFromL3/parallel_sort_l3.o
endif

ifeq (PM_FROM_L3,$(findstring PM_FROM_L3,$(CONFIGVARS))) # add L3 pm part
L3_OBJS  +=  ImportFromL3/pm_periodic_l3.o ImportFromL3/rfftwnd_mpi_threads.o ImportFromL3/transpose_mpi_threads.o
endif


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = P-Gadget3

OBJS  =  voronoi.o mpz_extension.o voronoi_density.o voronoi_force.o voronoi_2d.o voronoi_3d.o voronoi_exchange.o \
	 twopoint.o voronoi_meshrelax.o snIa_heating.o kspace_neutrinos.o \
	 cs_metals.o cs_sfr.o  cs_enrichment.o cs_findhotngbs.o cs_update_weights.o cs_ngbfind.o \
	 cosmic_rays_diffusion.o conduction.o fof.o subfind.o subfind_vars.o subfind_collective.o subfind_serial.o subfind_so.o subfind_cont.o \
	 subfind_distribute.o subfind_findlinkngb.o subfind_nearesttwo.o subfind_loctree.o subfind_alternative_collective.o subfind_reshuffle.o \
	 compare_partitions.o domain.o main.o greenf_diffusion.o driftfac.o growing_disk_potential.o \
	 subfind_potential.o subfind_density.o lineofsight.o kinfb_mhm.o sfr_mhm.o blackhole.o \
	 predict.o endrun.o global.o chemistry_noneq.o begrun.o run.o allvars.o\
	 timestep.o init.o restart.o io.o sfr_eff.o checksummed_sendrecv.o \
	 accel.o forcetree_update.o read_ic.o read_ic_cluster.o cooling.o ngb.o parallel_sort.o second_order.o \
	 system.o allocate.o density.o bubbles.o sizelimited_sendrecv.o rt_chem.o rt_bh_lum.o rt_sfr_lum.o rt_cooling.o \
	 gravtree.o kicks.o hydra.o rt_eddington.o rt_n.o rt_CGmethod.o rt_stars_lum.o rt_gas_lum.o  \
	 potential.o peano.o gravtree_forcetest.o parallel_sort_special.o \
	 pm_nonperiodic.o longrange.o mymalloc.o debug_md5/calc_checksum.o  debug_md5/Md5.o  \
	 cosmic_rays.o machfinder.o b_from_rot_a.o smooth_simple.o \
         phasespace.o phasespace_math.o helm_eos.o healpix.o compile_time_info.o \
         voronoi_ghost_search_alternative.o voronoi_ghost_search.o \
	 um_chemistry_noneq.o um_metal_cooling.o \
	 sinks.o \
	 cr_electrons.o \
	 scf.o scf_util.o \
	 mol_clouds.o

OBJS	+= $(GRAVITY_OBJS)

OBJS	+= $(L3_OBJS)

INCL   = cr_electrons.h cs_metals.h allvars.h proto.h forcetree.h cooling.h domain.h  cosmic_rays.h chemistry.h \
	 machfinder.h subfind.h dd.h fof.h voronoi.h helm_eos.h myqsort.h chemcool_consts.h kernel.h Makefile \

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS	+=  utilities.o integrate.o network_solver.o network.o 
INCL	+=  utilities.h integrate.h network_solver.h network.h
endif

ifeq (MODGRAV,$(findstring MODGRAV,$(CONFIGVARS)))
OBJS	+=  modgrav/modgrav_fieldsolve.o
INCL	+=  
endif

ifeq (WINDTUNNEL,$(findstring WINDTUNNEL,$(CONFIGVARS)))
OBJS	+=  windtunnel/tunnel.o
INCL	+=  
endif

ifeq (LT_METAL_COOLING,$(findstring LT_METAL_COOLING,$(CONFIGVARS)))
OBJS  += lt_wal_cooling.o lt_sn.o lt_io.o lt_utils.o lt_sfr.o ImfLib/libimf.o ImfLib/libimf_vars.o 
INCL  += lt_sn.h lt.h lt_error_codes.h lt_wal_cooling.h 
endif

ifeq (BP_REAL_CRs,$(findstring BP_REAL_CRs,$(CONFIGVARS))) # add bp cr part
OBJS += bp_cosmic_rays/bp_cosmic_rays.o
INCL += bp_cosmic_rays/bp_cosmic_rays.h
endif

ifeq (CHEMCOOL,$(findstring CHEMCOOL,$(CONFIGVARS)))
OBJS  += chemcool.o raytrace.o 
INCL  += chemcool_consts.h 
FOBJS = calc_photo.o calc_temp.o cheminmo.o compute_heating.o compute_md_kappa.o const_rates.o cool_func.o \
	cool_util.o coolinmo.o dvode.o evolve_abundances.o jac.o photoinit_lowZ.o \
	rate_eq_primordial.o spline.o validate_rates.o 
FINCL = cool.h fs_data.h mol_data.h non_eq.h shield_data.h Makefile 
else
FC    = $(CC)
endif

ifeq (EAGLE,$(findstring EAGLE,$(CONFIGVARS)))
OBJS  += eagle/eagle_cooling.o 
INCL  += eagle/eagle_cooling.h 
endif

ifeq (SAM_FEEDBACK,$(findstring SAM_FEEDBACK,$(CONFIGVARS)))
OBJS  += sam_feedback_vars.o sam_feedback_driver.o
INCL  += sam_feedback.h

ifeq (SAM_STARFORM,$(findstring SAM_STARFORM,$(CONFIGVARS)))
OBJS  += sam_starformation.o 
INCL  +=
endif

ifeq (SAM_SN2_FEEDBACK,$(findstring SAM_SN2_FEEDBACK,$(CONFIGVARS)))
ifeq (SN2_FEEDBACK_KINETIC,$(findstring SN2_FEEDBACK_KINETIC,$(CONFIGVARS)))
OBJS += sam_sn2_feedback_kinetic.o	
INCL +=
endif
endif

ifeq (SAM_AGN_FEEDBACK,$(findstring SAM_AGN_FEEDBACK,$(CONFIGVARS)))
ifeq (AGN_FEEDBACK_FIXED_RAD,$(findstring AGN_FEEDBACK_FIXED_RAD,$(CONFIGVARS)))
OBJS += sam_agn_feedback_fixed_rad.o	
INCL +=
endif
endif

ifeq (SAM_METALS,$(findstring SAM_METALS,$(CONFIGVARS)))
ifeq (METALS_FIXED_RAD,$(findstring METALS_FIXED_RAD,$(CONFIGVARS)))
	OBJS += sam_metals_fixed_rad.o
endif
endif

endif


CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL)

ifeq (VIP,$(findstring VIP,$(CONFIGVARS)))
FFLAGS = $(FOPTIONS)
else
FFLAGS = $(OPTIONS)
endif



FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 


LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW)

ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS))) 
LIBS   +=  -lpthread
endif

ifeq (VORONOI,$(findstring VORONOI,$(CONFIGVARS)))
LIBS   +=  $(GMP_LIBS) $(GMPLIB)
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c


$(FOBJS): $(FINCL)

complile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *~


