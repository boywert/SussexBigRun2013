EXEC   = L-Galaxies

OBJS   = ./code/main.o ./code/io_tree.o ./code/init.o ./code/cool_func.o \
     ./code/save.o ./code/save_galtree.o \
     ./code/mymalloc.o ./code/read_parameters.o \
	 ./code/peano.o ./code/allvars.o ./code/age.o ./code/update_type_two.o \
	 ./code/metals.o \
	 ./code/recipe_infall.o \
	 ./code/recipe_cooling.o \
	 ./code/recipe_starformation_and_feedback.o \
	 ./code/recipe_reincorporation.o \
	 ./code/recipe_mergers.o \
	 ./code/recipe_dust.o \
	 ./code/recipe_misc.o \
	 ./code/recipe_disrupt.o \
	 ./code/recipe_stripping.o \
	 ./code/scale_cosmology.o

INCL   = ./code/allvars.h  ./code/proto.h  Makefile

# Either include the default set of Makefile options, or define your own
#include Makefile_options
include My_Makefile_options

# Choose your system type (needs to match an entry in Makefile_compilers)
SYSTYPE = "Sussex"
include Makefile_compilers

LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas 

CFLAGS =   -g $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) My_Makefile_options Makefile_compilers

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	mkdir -p python
	${CC_MD} ${OPT} ${CFLAGS} -E ./code/allvars.h -o ./code/allvars.i
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/GALAXY_OUTPUT_2_TypeString.awk > ./awk/L-Galaxies_Types.txt
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/GALAXY_OUTPUT_2_DDL.awk > ./awk/L-Galaxies_DDL.sql	
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_struct.awk >  ./awk/idl/LGalaxy.pro
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/python/GALAXY_OUTPUT_2_python_struct.awk >  ./python/LGalaxyStruct.py
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_hists.awk > ./awk/idl/LGalaxy_plot.pro
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./awk/idl/LGalaxy_testfloats.pro
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/idl/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./awk/idl/LGalaxy_zerofloats.pro
	awk -f ./awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./awk/L-Galaxies.h
	awk -f ./awk/extractMOMAF_INPUT.awk ./code/allvars.i |awk -f ./awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./awk/L-Galaxies.h
