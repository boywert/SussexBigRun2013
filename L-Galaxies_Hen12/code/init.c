#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

/**@file init.c
 * @brief Sets up some unit conversion variables; converts SN and AGN feedback
 *        variables into internal units; and reads in input tables, including
 *        the desired output snapshots, the photometric and dust tables, the
 *        cooling functions and reionization tables.
 *
 *        <B>set_units()</B> - THIS IS FUNDAMENTAL TO UNDERSTAND THE UNITS IN
 *        THE CODE! sets ups some variables used to convert from physical to
 *        internal units (as UnitDensity_in_cgs). These are obtained from the
 *        three defined in input.par: UnitLength_in_cm, UnitMass_in_g and
 *        UnitVelocity_in_cm_per_s.
 *
 *        <B>read_output_snaps()</B> - reads in the list of output redshifts from
 *        file ./input/desired_output_redshifts.txt and converts them into snapsshots
 *        for the given cosmology.
 *
 *        <B>read_zlist()</B> - reads in 1/(z+1) from FileWithZList defined
 *        in ./input/input.par and creates a table with output
 *        redshift ZZ[] and ages Age[].
 *
 *        <B>read_recgas()</B> - Done if GASRECYCLE OFF - option not supported.
 *        Reads in Frac[] from gas_m62.dat.
 *
 *        <B>read_file_nrs()</B> - Done if SPECIFYFILENR OFF - the dark matter files
 *        to read can be defined in a file, instead of being read sequentially.
 *        These are defined in FileNrDir, in input.par, and read into
 *        ListInputFilrNr[].
 *
 *        <B>read_sfrz()</B> - Done if H2FORMATION OFF - option not supported.
 *        Reads in Rho[] and H2[][].
 *
 *        <B>read_reionization()</B> - Reads in Reion_z[] and Reion_Mc[] from
 *        ./input/Mc.txt. These are used if UPDATEREIONIZATION ON to get Okamoto(2008)
 *        fitting parameters (Mc), instead of Kravtsov(2004)  for the Gnedin (2000)
 *        reionization formalism.
 *
 *        <B>setup_spectrophotometric_model()</B> - Reads in the look up tables
 *        from a given Stellar Population Synthesis Model (**_Set_Phot_table_m**.dat).
 *        Default - Bruzual & Charlot 2003. Detailed description in
 *        add_to_luminosities() in recipe_misc.c. There are different files, each
 *        one corresponding to a different metallicity. On each file the tables have
 *        the Magnitudes for single bursts of 1\f$M_{\odot}\f$ for a range of ages
 *        and "inversely" k-corrected in case the code needs to compute observed
 *        frame magnitudes.
 *
 *
 *        <B>read_dust_tables()</B> - Reads in the dust extinction for the same bands
 *        read from the spectrophotometric tables both for the inter-galactic medium
 *        (**_Set_Ext_table.dat) and young birth clouds (**_Set_YSExt_table.dat).
 *        Detailed description at recipe_dust.c
 *
 *        <B>read_cooling_functions()</B> - calls the functions that read in the
 *        cooling functions defined in cool_func.c
 *
 *        In init.c, but called from other files, are the function to interpolate
 *        properties into different tables. They are all the same but interpolate
 *        into different properties (have different inputs):
 *
 *        <B>find_interpolation_point()</B> - Used by update_from_recycle() (not
 *        supported) to interpolate into the age table from the SSP look up tables.
 *
 *        <B>find_interpolated_lum()</B> - Used by add_to_luminosities() to
 *        interpolates into the age and metallicity in the SSP tables.
 *
 *        <B>find_interpolated_obslum()</B> - Same as find_interpolated_lum() but
 *        also interpolates in redshift. Not needed if the redshifts for the "inverted"
 *        kcorrections in the SSP tables correspond to the simulations redshifts.
 *
 *         <B>find_interpolate_reionization()</B> - Called from recipe_infall
 *        interpolates into the Reion_z[], Reion_Mc[] table, giving a value of Mc
 *        for a given redshift.
 *
 *        SuperNovae and AGN feedback parameters are converted into internal units:
 *
 *        \f$ AgnEfficiency = \frac{UnitMass_{\rm{g}}}{1.58e^{-26}UnitTime_{\rm{s}}}\f$
 *
 *        \f$ EnergySNcode = \frac{EnergySN}{UnitEnergy_{\rm{cgs}}} h; \f$
          \f$ EtaSNcode = EtaSN \frac{UnitMass_{\rm{g}}}{M_\odot h}. \f$

 *          */

//Needs to be moved to proto.h
void read_reionization(void);

/**@brief controlling recipe for init.c, calls functions to read in tables and
 *        defines some SN and AGN feedback parameters.
 *
 *        Calls set_units(), read_output_snaps(), read_zlist(), read_recgas(),
 *        read_file_nrs(), read_sfrz(), read_reionization(),
 *        setup_spectrophotometric_model(), read_dust_tables() and
 *        read_cooling_functions().
 *
 *        Converts EnergySN (->EnergySNcode), EtaSN (->EtaSNcode) and AgnEfficiency
 *        into internal units.
 *        */
void init(void)
{
  int i;
  struct rlimit rlim;

  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

// define log of global metalicities
  for(i = 0; i < IZZ; i++)
    {
      log10zz[i] = log10(def_zz[i]);
    }

  set_units();
#ifdef GALAXYTREE
  ScaleFactor = pow(2, Hashbits) / BoxSize;
#endif

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
  EtaSNcode = EtaSN * (UnitMass_in_g / SOLAR_MASS) / Hubble_h;
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("\nEnergySNcode*EtaSNcode= %g\n", EnergySNcode * EtaSNcode);
#else
  printf("\nEnergySNcode*EtaSNcode= %g\n", EnergySNcode * EtaSNcode);
#endif

  //reads in the redshifts for the used Cosmology
  read_zlist();
  //reads in the redshifts for Original Cosmology
  read_zlist_original_cosm();

  //reads in the desired output snapshots
   read_output_snaps();


#ifdef GASRECYCLE
  /* read in the table of gas recycling */
  read_recgas();
#endif

#ifdef SPECIFYFILENR
  /* read in the number of the files to be processed */
  read_file_nrs();
#endif

#ifdef H2FORMATION
  read_sfrz();
#endif

  if(ReionizationOn == 2) 
    read_reionization();


  //Values of a for the beginning and end of reionization
  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

//CREATE ARRAYS OF SFH TIME STRUCTURE:
#ifdef  STAR_FORMATION_HISTORY
  create_sfh_bins();
#endif

//read in photometric tables
  setup_spectrophotometric_model();

  read_cooling_functions();

//READ IN THE YIELD TABLES, AND FORM NORMALISED YIELD ARRAYS:
#ifdef YIELDS
  read_yield_tables();
  integrate_yields();
#endif

  /* TODO - this conversion is done again in recipe_cooling.c - */


  /* Msun/year into g/s and then into internal units
   * of 1e10Msun/h / Mpc/Km/s/h */
  AgnEfficiency /= (UnitMass_in_g / UnitTime_in_s * 1.58e-26);


}



/*@brief Reads in 1/(z+1) from FileWithZList defined
 *       in ./input/input.par for the list of output snapshots.
 *       Creates a table with redshift ZZ[] and ages Age[].*/
void read_zlist(void)
{
	int i;
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList);

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      exit(1);
    }

  Zlistlen = 0;
  do
    {
      if(fscanf(fd, " %lg ", &AA[Zlistlen]) == 1)
      	Zlistlen++;
      else
      	break;
    }
  while(Zlistlen < MAXSNAPS);



  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
     {
       //convert AA[] into redshift - ZZ[]
       ZZ[i] = 1 / AA[i] - 1;
       //table with time in internal units (Mpc/Km/s/h)
       if(ZZ[i]>=0.0)
         Age[i] = time_to_present(ZZ[i]);
       else
      	 Age[i] = 0.0;
    	 //break;
     }

#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}


void read_zlist_new(void)
{
  int i, dummy_snap;
  double dummy_a;
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList);
  if(!(fd = fopen(fname, "r")))
    {
  	char sbuf[1000];
  	sprintf(sbuf, "can't open file `%s'\n", fname);
  	terminate(sbuf);
    }

  Zlistlen = 0;
  do
    {
      if(fscanf(fd, "%d %lg %lg", &dummy_snap, &ZZ[Zlistlen], &dummy_a) == 3)
    	  Zlistlen++;
      else
	break;

    }
  while(Zlistlen < MAXSNAPS);
  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
  {
  	//convert redshift - ZZ[] into AA[]
  	AA[i] = 1/(ZZ[i] + 1);
  	//printf("z[%d]=%f\n",i,ZZ[i]);
  	//table with time in internal units (Mpc/Km/s/h)
  	if(ZZ[i]>=0.0)
  		Age[i] = time_to_present(ZZ[i]);
  	else
  		Age[i] = 0.0;
  		//break;
  }


#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}

void read_zlist_original_cosm(void)
{
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList_OriginalCosm);
  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      exit(1);
    }

  Zlistlen = 0;
  do
    {
      if(fscanf(fd, " %lg ", &AA_OriginalCosm[Zlistlen]) == 1)
        {Zlistlen++;
	}
      else
        break;
    }
  while(Zlistlen < MAXSNAPS);

  fclose(fd);

#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
}


/**@brief Reads in the list of output snapshots from
 *        file /input/desired_output_snaps.txt*/
void read_output_snaps(void)
{
  int i, j;

#ifndef GALAXYTREE
  char buf[1000];
  FILE *fd;

  LastSnapShotNr=0;

  sprintf(buf, "%s", FileWithOutputRedshifts);

  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  for(i = 0; i < NOUT; i++)
    {
      if(fscanf(fd, " %f ", &ListOutputRedshifts[i]) != 1)
        {
    	  char sbuf[1000];
    	  sprintf(sbuf, "I/O error in file '%s'\n", buf);
    	  terminate(sbuf);
        }

      //find the snapshot corresponding to the desired output redshift in ListOutputRedshifts[]
      for(j = 0; j < MAXSNAPS; j++)
    	  if(ListOutputRedshifts[i]>=ZZ[j])
    	    {
    		  if((ZZ[j-1]-ListOutputRedshifts[i])<(ListOutputRedshifts[i]-ZZ[j]) || ZZ[j]< 0.0)
    	        ListOutputSnaps[i]=j-1;
    		  else
    		    ListOutputSnaps[i]=j;
    		  //printf("outz=%f zz[%d]=%f snap=%d\n",ListOutputRedshifts[i], j-1, ZZ[j-1]);
    		  break;
    	    }

      //define the LastSnapShotNr as the highest snapshot need to be computed
      if(LastSnapShotNr<ListOutputSnaps[i])
    	LastSnapShotNr=ListOutputSnaps[i];
    }
  fclose(fd);




#else
  for(i = 0; i < NOUT; i++)
    ListOutputSnaps[i] = i;
  LastSnapShotNr=LastDarkMatterSnapShot;
#endif
}


/**@brief Reads in the look up tables from Stellar Population Synthesis Models.
 *
 * Reads in the look up tables from a given Stellar Population Synthesis Model.
 * The default model is Bruzual & Charlot 2003. There are different files, each
 * one corresponding to a different metallicity. On each file the tables have
 * the Magnitudes for single bursts of 1\f$M_{\odot}\f$ for a range of ages
 * and "inversely" k-corrected in case the code needs to compute observed
 * frame magnitudes. 6 different metallicities are available at ./PhotTables/
 * with the file names **_Set_Phot_table_m22.dat/m32/m42/m52/m62/m72
 * corresponding to Z=0.0001/z=0.0004/z=0.004/Z=0.008/Z=0.02/Z=0.05/Z=0.1. Two
 * different prefixes are available DataRelease or SDSS corresponding respectively
 * to BVRIK or ugriz bands (both with all the different metallicitites.)
 *
 * On each file, the three first numbers correspond to: the number of snapshots,
 * the number of bands and the number of age bins. Then the age grid is listed.
 * After that for each snapshot, for each age the value corresponding to a burst
 * inversely k-corrected to that redshift, with that age (and metallicity) is listed.
 *
 * The basic structure is number of columns  corresponding to number of mags,
 * repeated over the number of ages on the initial listed grid, multiplied by the
 * number of snapshots, with the snapshot number listed in between.
 *
 * agTableZz[Mag][mettallicity][Snapshot][Age]. */
void setup_spectrophotometric_model(void)
{
  FILE *fa, *fb;
  int i, k, p, band;
  char buf[100], FilterName[100], dummy[100];
  int nsnaps, n_age_bins;
#ifdef MRII
  char Simulation[100] = { "MRII" };
#else
  char Simulation[100] = { "MR" };
#endif
#ifdef M05
  char *filenames[] = { "m0001", "m001", "m002", "m004" };
#else
  char *filenames[] = { "m00001", "m00004", "m0004", "m0008", "m002", "m005" };
#endif

  /*Loop over the different files corresponding to different metallicities */
  for(k = 0; k < IZZ; k++)
    {
      sprintf(buf, "%s", FileWithFilterNames);
      if((fa = fopen(buf, "r")) == NULL)
	printf("\n**%s not found on line %d of init.c**\n", buf, __LINE__);

      FilterLambda[NMAG] = 0.55;	//to use by the dust model, the wavelength of the V-band
      //There is a different file for each band
      for(band = 0; band < NMAG; band++)
	{
	  fscanf(fa, "%s", FilterName);


	  sprintf(buf, "%s/%s_Phot_Table_%s_Mag%s_%s.dat", PhotDir, PhotPrefix, Simulation, FilterName,
		  filenames[k]);
	  if(!(fb = fopen(buf, "r")))
	    {
	      char sbuf[1000];
	      sprintf(sbuf, "file `%s' not found.\n", buf);
	      terminate(sbuf);
	    }
#ifdef PARALLEL
	  if(ThisTask == 0)
	    printf("reading file %s \n", buf);
#else
	  printf("reading file %s \n", buf);
#endif
	  fscanf(fb, "%s %f %d %d", dummy, &FilterLambda[band], &nsnaps, &n_age_bins);
	  /* check that the numbers on top of the file correspond to NMAG and ZL_LEN */
	  if(nsnaps != IZ)
	    {
	      terminate("nsnaps not equal to  Snaplistlen");
	    }
	  if(n_age_bins != ZL_LEN)
	    {
	      terminate("n_age_bins not equal to ZL_LEN");
	    }

	  // read ages of SSPs
	  for(p = 0; p < ZL_LEN; p++)
	    {
	   	   fscanf(fb, " %e ", &AgeTab[p]);
	   	   if(AgeTab[p] > 0.0)	// avoid taking a log of 0 ...
		     {
			   /* converts AgeTab from years to log10(internal time units 1e12 Yrs/h)
			   *UnitTime_in_Megayears is ~1e6, but not quite. look at set_units() at init.c*/
			   AgeTab[p] = AgeTab[p] / 1.0e6 / UnitTime_in_Megayears * Hubble_h;
			   AgeTab[p] = log10(AgeTab[p]);
			 }
	    }

	  // read mags at each output redshift
	  for(p = 0; p < IZ; p++)
	    {
	      fscanf(fb, " %f ", &RedshiftTab[p]);
	      //for each age
	      for(i = 0; i < ZL_LEN; i++)
			{
			  fscanf(fb, "%e", &MagTableZz[band][k][p][i]);
			  MagTableZz[band][k][p][i] = pow(10., -MagTableZz[band][k][p][i] / 2.5);
			}		//end loop on age
		}		//end loop on redshift (everything done for current band)

	  fclose(fb);
	}			//end loop on bands

      fclose(fa);
    }				//end loop on metallicity

  init_jump_index();
}


#define NJUMPTAB 1000
int jumptab[NJUMPTAB];
double jumpfac;

void init_jump_index(void)
{
  double age;
  int i, idx;

  jumpfac = NJUMPTAB / (AgeTab[ZL_LEN - 1] - AgeTab[1]);

  for(i = 0; i < NJUMPTAB; i++)
    {
      age = AgeTab[1] + i / jumpfac;
      idx = 1;
      while(AgeTab[idx + 1] < age)
	idx++;
      jumptab[i] = idx;
    }
}


int get_jump_index(double age)
{
  return jumptab[(int) ((age - AgeTab[1]) * jumpfac)];
}

/**@brief Used by update_from_recycle() (notsupported) to interpolate
 *        into the age table from the SSP look up tables.*/
void find_interpolation_point(double timenow, double timetarget, int *index, double *f1, double *f2)
{
  int idx;
  double age, frac;

  age = timenow - timetarget;

  if(age > 0)
    {
      age = log10(age);

      if(age > AgeTab[ZL_LEN - 1])	/* beyond table, take latest entry */
	{
	  *index = ZL_LEN - 2;
	  *f1 = 0;
	  *f2 = 1;
	}
      else if(age < AgeTab[1])	/* age younger than 1st enty, take 1st entry */
	{
	  *index = 0;
	  *f1 = 0;
	  *f2 = 1;
	}
      else
	{
	  /*
	     idx = 1;
	   */
	  idx = get_jump_index(age);
	  while(AgeTab[idx + 1] < age)
	    idx++;
	  *index = idx;
	  frac = (age - AgeTab[idx]) / (AgeTab[idx + 1] - AgeTab[idx]);
	  *f1 = 1 - frac;
	  *f2 = frac;
	}
    }
  else				/* this lies in the past */
    {
      *index = 0;
      *f1 = 0;
      *f2 = 0;
    }

}

/**@brief Used by add_to_luminosities() to interpolates into
 *        the age and metallicity in the SSP tables.*/
void find_interpolated_lum(double timenow, double timetarget, double metallicity, int *metindex,
			   int *tabindex, double *f1, double *f2, double *fmet1, double *fmet2)
{
  int k, i, idx;
  double age, frac;
  double ft1, ft2, fm1, fm2;

  age = timenow - timetarget;

  if(age > 0)
    {
      age = log10(age);

      if(age > AgeTab[ZL_LEN - 1])	/* beyond table, take latest entry */
	{
	  k = ZL_LEN - 2;
	  ft1 = 0;
	  ft2 = 1;
	}
      else if(age < AgeTab[1])	/* age younger than 1st enty, take 1st entry */
	{
	  k = 0;
	  ft1 = 0;
	  ft2 = 1;
	}
      else
	{
	  /* 
	     idx = 1; 
	   */
	  idx = get_jump_index(age);
	  while(AgeTab[idx + 1] < age)
	    idx++;
	  k = idx;
	  frac = (age - AgeTab[idx]) / (AgeTab[idx + 1] - AgeTab[idx]);
	  ft1 = 1 - frac;
	  ft2 = frac;
	}
    }
  else				/* this lies in the past */
    {
      k = 0;
      ft1 = 0;
      ft2 = 0;
    }

  /* Now interpolate also for the metallicity */
  metallicity = log10(metallicity);

  if(metallicity > log10zz[IZZ - 1])	/* beyond table, take latest entry */
    {
      i = IZZ - 2;
      fm1 = 0;
      fm2 = 1;
    }
  else if(metallicity < log10zz[0])	/* mettallicity smaller 1st enty, take 1st entry */
    {
      i = 0;
      fm1 = 1;
      fm2 = 0;
    }
  else
    {
      idx = 0;
      while(log10zz[idx + 1] < metallicity)
	idx++;
      i = idx;
      frac = (metallicity - log10zz[idx]) / (log10zz[idx + 1] - log10zz[idx]);
      fm1 = 1 - frac;
      fm2 = frac;
    }


  *metindex = i;
  *tabindex = k;

  *f1 = ft1;
  *f2 = ft2;
  *fmet1 = fm1;
  *fmet2 = fm2;
}


/**@brief Same as find_interpolated_lum() but also interpolates in redshift;
 *        Not needed if the redshifts for the "inverted" kcorrections in the
 *        SSP tables correspond to the simulation's redshifts.*/
void find_interpolated_obslum(double timenow, double timetarget, double metallicity, double redshift,
			      int *metindex, int *tabindex, int *redindex, double *f1, double *f2,
			      double *fmet1, double *fmet2, double *fred1, double *fred2)
{
  int k, i, j, idx;
  double age, frac;
  double ft1, ft2, fm1, fm2, fz1, fz2;

  age = timenow - timetarget;

  if(age > 0)
    {
      age = log10(age);

      if(age > AgeTab[ZL_LEN - 1])	/* beyond table, take latest entry */
	{
	  k = ZL_LEN - 2;
	  ft1 = 0;
	  ft2 = 1;
	}
      else if(age < AgeTab[1])	/* age younger than 1st enty, take 1st entry */
	{
	  k = 0;
	  ft1 = 0;
	  ft2 = 1;
	}
      else
	{
	  /* 
	     idx = 1; 
	   */
	  idx = get_jump_index(age);
	  while(AgeTab[idx + 1] < age)
	    idx++;
	  k = idx;
	  frac = (age - AgeTab[idx]) / (AgeTab[idx + 1] - AgeTab[idx]);
	  ft1 = 1 - frac;
	  ft2 = frac;
	}
    }
  else				/* this lies in the past */
    {
      k = 0;
      ft1 = 0;
      ft2 = 0;
    }

  /* Now interpolate also for the metallicity */
  metallicity = log10(metallicity);

  if(metallicity > log10zz[IZZ - 1])	/* beyond table, take latest entry */
    {
      i = IZZ - 2;
      fm1 = 0;
      fm2 = 1;
    }
  else if(metallicity < log10zz[0])	/* age younger than 1st enty, take 1st entry */
    {
      i = 0;
      fm1 = 1;
      fm2 = 0;
    }
  else
    {
      idx = 0;
      while(log10zz[idx + 1] < metallicity)
	idx++;
      i = idx;
      frac = (metallicity - log10zz[idx]) / (log10zz[idx + 1] - log10zz[idx]);
      fm1 = 1 - frac;
      fm2 = frac;
    }

  /* Now interpolate also for the redshift of the observations */

  if(redshift > RedshiftTab[IZ - 1])	/* beyond table, take latest entry */
    {
      j = IZ - 2;
      fz1 = 0;
      fz2 = 1;
    }
  else if(redshift < RedshiftTab[0])	/* age younger than 1st enty, take 1st entry */
    {
      j = 0;
      fz1 = 1;
      fz2 = 0;
    }
  else
    {
      idx = 0;
      while(RedshiftTab[idx + 1] < redshift)
	idx++;
      j = idx;
      frac = (redshift - RedshiftTab[idx]) / (RedshiftTab[idx + 1] - RedshiftTab[idx]);
      fz1 = 1 - frac;
      fz2 = frac;
    }


  *metindex = i;
  *tabindex = k;
  *redindex = j;

  *f1 = ft1;
  *f2 = ft2;
  *fmet1 = fm1;
  *fmet2 = fm2;
  *fred1 = fz1;
  *fred2 = fz2;
}


/**@brief NEVER USED*/
int find_index(float *xx, int n, float x)
{
//TODO - never used
  int j, jl, ju, jm;

  jl = 0;
  ju = n + 1;

  while((ju - jl) > 1)
    {

      jm = (ju + jl) / 2.;

      if((xx[n - 1] > xx[0]) && (x > xx[jm - 1]))
	jl = jm;

      else
	ju = jm;

    }

  j = jl;
  if(j == 0)
    j = j + 1;
  return j;
}

/**@brief Sets up some variables used to convert from physical to internal
 *        units (as UnitDensity_in_cgs); These are obtained from the three
 *        defined in input.par: UnitLength_in_cm (cm to Mpc), UnitMass_in_g
 *        (g to 1e10Msun) and UnitVelocity_in_cm_per_s (cm/s to km/s).
 *
 *       As defined in input.par, \f$\rm{UnitLength}_{\rm{cm}}=
 *       3.08568\times 10^{24}\rm{cm}\f$, converts from cm into Mpc and
 *       \f$\rm{UnitVelocity}_{\rm{cm/s}}=10000\rm{cm/s}\f$, converts from
 *       cm/s to Km/s (cm to km). In set_units() \f$\rm{UnitTime}_{\rm{s}}\f$
 *       is derived from these two quantities:
 *
 *       \f$\frac{\rm{UnitLength}_{\rm{cm}}}{\rm{UnitVelocity}_{\rm{cm/s}}}
 *       =3.08568\times 10^{19}\rm{Mpc}~\rm{Km}^{-1}\rm{s}^{-1}\f$,
 *
 *       through the code \f$t_{\rm{dyn}}\f$ has internal units and its never
 *       converted (note that \f$t_{\rm{dyn}}\f$ has an h factor, as the code internal
 *       units which despite not being included is \f$\rm{UnitTime}_{\rm{s}}\f$ is
 *       included in the output of time_to_present() - so it is consistent).
 *
 *       \f$ \rm{UnitDensity}_{\rm{cgs}} =
 *       \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}}^3}=6.769898\times 10^{-31}\f$,
 *       converts density in \f$\rm{g}~\rm{cm}^{-3}\f$ into internal units
 *       \f$(10^{10}M_{\odot}\rm{Mpc}^{-3})\f$
 *
 *        \f$ \rm{UnitPressure}_{\rm{cgs}} =
 *        \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}} \times \rm{UnitTime}_{\rm{s}}^2}
 *        =6.769898\times 10^{-21}\f$, converts pressure in
 *        \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-2}\f$ into internal units
 *        \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s) \f$
 *
 *       \f$ \rm{UnitCoolingRate}_{\rm{cgs}} =
 *       \frac{\rm{UnitPressure}_{\rm{cgs}}}{\rm{UnitTime}_{\rm{s}}}=2.193973\times 10^{-40}\f$,
 *       converts the cooling rate in \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-3}\f$ into
 *       internal units \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s)^{-3}) \f$
 *
 *      \f$ \rm{UnitEnergy}_{\rm{cgs}} =
 *       \frac{\rm{UnitMass}_{\rm{g}} \times \rm{UnitLength}_{\rm{cm}}^2}{\rm{UnitTime}_{\rm{s}}^2}
 *       =1.989000\times 10^{53}\f$, converts energy in
 *       \f$\rm{g}~\rm{cm}^2\rm{s}^{-2}\f$ into internal units
 *       \f$(10^{10}M_{\odot}~\rm{Mpc}^{2}(Mpc/Mk/s)^{-2})\f$
 *
 *       \f$ \rm{Hubble} = \rm{HUBBLE} \times \rm{UnitTime}_{\rm{s}}=100.0001\f$, where
 *       \f$\rm{HUBBLE}=3.2407789\times 10^{-18} h~\rm{s}^{-1}\f$, is the hubble
 *       constante in \f$(h~\rm{Km}~\rm{s}^{-1}\rm{Mpc}^{-1})\f$.
 *
 *       Then define a constant:
 *       \f$ \rm{RhoCrit} = \frac{3\times \rm{Hubble}^2}{8\times \rm{M}\_\rm{PI} \times G}=
 *       27.75505~h^2~10^{10}M_{\odot}\rm{Mpc}^{-3}\f$,
 *
 *       */

void set_units(void)
{


	// SEC_PER_MEGAYEAR   3.155e13
	// SEC_PER_YEAR       3.155e7

  //UnitLength_in_cm & UnitVelocity_in_cm_per_s; defined at input.par
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  UnitTime_in_years = 1e6*UnitTime_in_Megayears;


  //gravity in internal units
  G = GRAVITY / pow3(UnitLength_in_cm) * UnitMass_in_g * pow2(UnitTime_in_s);//43.00708

  //converts g.cm^-3 into internal units (1e10Msun Mpc^-3)
  UnitDensity_in_cgs = UnitMass_in_g / pow3(UnitLength_in_cm);//6.769898e-31

  //converts g.cm^-1s^-2 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-2) \f$
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow2(UnitTime_in_s);//6.769898e-21

  //converts g.cm^-1.s^-3 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-3)
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;//2.193973e-40

  //converts g.cm^2.s^-2 into internal units (10^10Msun.Mpc^2(Mpc/Km/s)^-2)
  UnitEnergy_in_cgs = UnitMass_in_g * pow2(UnitLength_in_cm) / pow2(UnitTime_in_s);//1.989000e+53


  /* converts the Hubble constant from h.s^-1 into h.Km.s^-1.Mpc-1 */
  // Would make much more sense to define this including the h100 factor.
  Hubble = HUBBLE * UnitTime_in_s;//100.000

  /* compute a quantitity */
//TODO - never used
  RhoCrit = 3 * Hubble * Hubble / (8 * M_PI * G);//27.75505 (h^2.10^10Msun.Mpc^-3)

}

#ifdef GASRECYCLE
void read_recgas(void)
{
  char buf[1000];
  FILE *fd;
  int p, i, k, j, l, dummy;
  float dumb;
  if(!(fd = fopen("./gas_m62.dat", "r")))
    terminate("file not found.\n");
    

  for(p = 0; p < ZL_LEN; p++)
    {
      fscanf(fd, "%f %f", &dumb, &Frac[p]);
      /*     printf("Frac reading: dumb %f , frac%f\n",dumb,Frac[p]); */
    }
  fclose(fd);

}
#endif


#ifdef SPECIFYFILENR

void read_file_nrs(void)
{
  int i;
  char buf[1000];
  FILE *fd;
  sprintf(buf, "%s", FileNrDir);
  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  for(i = 0; i < 111; i++) //only 111files in ../input/filenrdir.txt are read in
    {
      if(fscanf(fd, " %d ", &ListInputFilrNr[i]) != 1)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "I/O error in file '%s'\n", buf);
	  terminate(sbuf);
	}
    }
  fclose(fd);
  
}


#endif

#ifdef H2FORMATION

void read_sfrz(void)
{
  char buf[1000];
  FILE *fd;
  int p, i, k, j, l, dummy;
  float dumb;
  if(!(fd = fopen("./SFR_Z.dat","r")))
    terminate("file not found.\n");
    
  printf("rho_len %d, zlen %d \n", RHO_LEN, Z_LEN);
  for(p = 0; p < RHO_LEN; p++)
    {
      fscanf(fd, "%f", &Rho[p]);
      for(i = 0; i < Z_LEN; i++)
	fscanf(fd, "%f", &H2[p][i]);

    }

  fclose(fd);

}

void find_interpolate_h2(double metalicity, double rho, int *tabindex, int *metindex, double *f1, double *f2,double *fmet1, double *fmet2)
{
  double frac;
  int idx;
  float zz[] = { -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1. };

  if (rho < Rho[0]) /* rho less than the smalles rho in the table, take the 1st entry*/
    {
      *tabindex = 0;
      *f1 = 1.;
      *f2 = 0.;
    }
  else if(rho > Rho[RHO_LEN - 1])	/* rho great than the largest rho in the table, take the last entry */
    {
      *tabindex = RHO_LEN - 2;
      *f1 = 0.;
      *f2 = 1.;
    }
  else
    {
      idx = 0;
      while(Rho[idx + 1] < rho)
	idx++;

      frac = (rho - Rho[idx]) / (Rho[idx + 1] - Rho[idx]);
      *f1 = 1 - frac;
      *f2 = frac;
      *tabindex = idx;
    }

  if(metalicity < zz[0])	/* rho less than the smalles rho in the table, take the 1st entry */
    {
      *metindex = 0;
      *fmet1 = 1.;
      *fmet2 = 0.;
    }
  else if(metalicity > zz[Z_LEN - 1])	/* rho great than the largest rho in the table, take the last entry */
    {
      *metindex = RHO_LEN - 2;
      *fmet1 = 0.;
      *fmet2 = 1.;
    }
  else
    {
      idx = 0;
      while(zz[idx + 1] < metalicity)
	idx++;

      frac = (-metalicity + zz[idx]) / (-zz[idx + 1] + zz[idx]);
      *fmet1 = 1 - frac;
      *fmet2 = frac;
      *metindex = idx;
    }

}
#endif


void read_reionization(void)
{
  FILE *fd;
  int p;
  float dumb;

  if(!(fd = fopen(McFile, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", McFile);
      terminate(sbuf);
    }

  for(p = 0; p < 45; p++)
    {
      fscanf(fd, "%f", &dumb);
      fscanf(fd, "%f", &Reion_z[p]);
      fscanf(fd, "%f", &Reion_Mc[p]);

    }
  Reion_z[45] = Reion_z[44];
  Reion_Mc[45] = Reion_Mc[44];


  fclose(fd);


}

void find_interpolate_reionization(double zcurr, int *tabindex, double *f1, double *f2)
{
  double frac;
  int idx;

  if(zcurr > Reion_z[0])	/* redshift is higher than the largest z in the table, take the 1st entry */
    {
      *tabindex = 0;
      *f1 = 1.;
      *f2 = 0.;
    }
  else if(zcurr <= Reion_z[44])	/* redshift smaller than the smallest rho in the table, take the last entry */
    {
      *tabindex = 44;
      *f1 = 1.;
      *f2 = 0.;
    }
  else
    {
      idx = 0;
      while(Reion_z[idx + 1] > zcurr)
	idx++;

      frac = (-zcurr + Reion_z[idx]) / (Reion_z[idx] - Reion_z[idx + 1]);
      *f1 = 1 - frac;
      *f2 = frac;
      *tabindex = idx;
    }

}
