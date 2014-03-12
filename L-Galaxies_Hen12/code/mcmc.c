#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-EPS)
#define EPS 1.2e-7


/**@file mcmc.c
 * @brief This is the main mcmc file, it reads in observational tests,
 *        starts a chain, proposes new sets of parameters, calls the
 *        SAM, gets the likelihood for galaxy properties obtained with
 *        the current parameters, compares the likelihood with the
 *        previous run and decides on whether or not to accept the
 *        proposed parameters.
 *
 * To run you need to chose input_mcmc_wmp7.par and select the
 * desired sample of halos. Also choose the correct desired_output_snap
 * according to the redshift at which you want to run the constraints.
 * Turn on MCMC on the Makefile.
 *
 * Senna is the controlling routine. First is calls read_observations,
 * to read in all the observational data sets into the struct MCMC_Obs.
 * Then, an initial SAM is run, to obtain a first likelihood value.
 * This will be the starting point of the chain, given by the initial
 * values of p[] (this should be some non 0 likelihood region to avoid
 * a long burn in). The likelihood for each set of parameters is computed
 * in mcmc_likelihood.c (get_likelihood()) at the end of each SAM run.
 *
 * After this first step, a chain is started with mcmc() being called
 * for each step to control all the processes. Namely, propose_new_parameters(),
 * SAM() and then at each step the likelihood for the given set of
 * parameters compared with previous and accepted or not. Two options
 * available according to the value of MCMCMode in input.par. If MCMCMode
 * =0, normal MCMC is done. If MCMCMode = 1, the new set of parameters
 * is only accepted if the likelihood is higher than for the previous.
 * This makes the chain going up in likelihood very quickly and its
 * useful to find a region of non 0 likelihood when the parameter is
 * still unknown.
 *
 * The name of the parameters sampled can be seen in the beginning of
 * function SAM() where the values of the internal variables of the code
 * containing the parameters are changed according to the new set of
 * parameters proposed. At each step p[] will contain the previously
 * accepted set of parameters and prop[] will contain the newly proposed
 * set of parameters (given by propose_new_parameters()). If desired,
 * less parameters then the default number, given by Npar can be sampled.
 * This can be done in propose_new_parameters() by never assigning a new
 * value to a given parameter. */
#ifdef MCMC
void Senna()
{
  int i, j, snap, chainweight = 1, IndividualAcceptRate = 0, TotAcceptRate = 0;
  FILE *fa;
  char buf[1000];
  char DirNum;

  printf("\n\n\n");
  printf("**********************************************************\n");
  printf("*                                                        *\n");
  printf("*                   Starting Senna                       *\n");
  printf("*                                                        *\n");
  printf("*             MCMC parameter estimation                  *\n");
  printf("*  Applied to a Semi-Analytic Model of Galaxy Formation  *\n");
  printf("*                                                        *\n");
  printf("**********************************************************\n\n");

  initialize_mcmc_par ();

  time_t initial, final;

  //This file will have the values for the parameters
  //It will be the output fom the mcmc sampling used in the analysis
  getcwd(buf, sizeof(buf));
  DirNum=buf[strlen(buf)-1];
  sprintf(buf, "%ssenna_g%c_%d.txt", OutputDir, DirNum, ThisTask);
  if((fa = fopen(buf, "w")) == NULL)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}


  //Read observational constraints (SMF, LFs, BHBM relation, etc)
  read_observations();

  //runs the SAM with the previous set of parameters and returns a likelyhood
  printf("\nTask %d running initial SAM\n\n",ThisTask);
  time(&initial);
  lhood1 = SAM(MCMCTreeSampleFile);
  time(&final);
  printf("Initial SAM took %lds\n", final - initial);

  MCMCseed = -(ThisTask * 100 + 25);
  for(i=1;i<ChainLength+1;i++)
	{
	  CurrentMCMCStep=i;
	  time(&initial);
	  printf("\n\n\nMCMC %d STARTED\n\n\n", i);
	  //returns either 0 or 1 and calculates lhood1
	  IndividualAcceptRate = mcmc(1);
	  TotAcceptRate = TotAcceptRate + IndividualAcceptRate;

	  if(lhood1<1e-99)
	  	lhood1 = 1e-99;
	  fprintf(fa," %d %0.10g",chainweight, -(log10(lhood1)));

	  //print parameters into output file
	  for(j=0;j<Npar;++j)
	  {
	  	if(MCMC_PAR[j].Sampling_Switch==1)
	  	{
	  		if(strcmp(MCMC_PAR[j].Type,"Physical")==0)
	  		{
	  			if(Sample_Physical_Parameters==1)
	  			{
	  				if(Time_Dependant_PhysPar==1)
	  					for(snap=0;snap<NOUT;snap++)
	  						fprintf(fa," %0.10g", log10(MCMC_PAR[j].Value[snap]));
	  				else
	  					fprintf(fa," %0.10g", log10(MCMC_PAR[j].Value[0]));
	  			}
	  		}
	  		else if (strcmp(MCMC_PAR[j].Type,"Cosmological")==0)
	  		{
	  			if(Sample_Cosmological_Parameters==1)
	  				fprintf(fa, " %0.10g", log10(MCMC_PAR[j].Value[0]));
	  		}
	  	}
	  }
	  fprintf(fa,"\n");
      fflush(fa);
	  time(&final);
	  printf("Task %d chain %d took %lds", ThisTask, i, final - initial);
	}
  fclose(fa);

  /* For each set of parameters writes the semi-analytic
   * predictions used for the likelihood calculation*/
  write_comparison_to_obs();

  myfree(MCMC_Obs);


  printf("\nAcceptance rate of this chain=%f%%\n", ((float) TotAcceptRate / ChainLength) * 100);
  printf("\n\nMCMC OVER\n\n");

}


////////
//MCMC//
////////


int  mcmc (int Nsteps)
{
  double lhood2, qratio, AcceptanceProbability, ran;
  int AcceptanceCounter = 0, i, AcceptanceLogic, snap;

  for(i = 0; i < Nsteps; i++)
	{
	  //get a new set of parameters and return qratio - the prior
	  qratio = propose_new_parameters();

	  //runs the SAM with the new parameters and gives the likelyhood for them
	  lhood2=SAM(MCMCTreeSampleFile);

	  printf("LIKELY1=%0.5e\nLIKELY2=%0.5e\n",lhood1,lhood2);

	  if(isnan(lhood1))
		  lhood1=0.0;
	  if(isnan(lhood2))
		  lhood2=0.0;

	  /* By default qratio = 1, meaning we assume a flat prior. Therefore,
	   * the acceptance probability is just given by the ratio of the likelihoods
	   * from two steps. */
	  if(1.0 > (qratio * (lhood2 / lhood1)))
		AcceptanceProbability = qratio * (lhood2 / lhood1);
	  else
		AcceptanceProbability = 1.0;

	  //reset seed as well
	  ran = ran3(&MCMCseed);

	  //accepts the proposed parameters with probability=acceptance rate
	  if(MCMCMode == 0)
		AcceptanceLogic = (ran < AcceptanceProbability);
	  //only accepts new parameters if likelihood increases
	  if(MCMCMode == 1)
		AcceptanceLogic = (lhood2 > lhood1);

	  //AcceptanceLogic=1;
	  if(AcceptanceLogic)
		{
		  for(i=0;i<Npar;++i)
		    for(snap=0;snap<NOUT;snap++)
		    	MCMC_PAR[i].Value[snap] = MCMC_PAR[i].PropValue[snap];

		  lhood1=lhood2;
		  AcceptanceCounter++;
		  printf("\n******************************************************\n");
		  printf("Accepted!!!\n");
		  for(i=0;i<Npar;++i)
		  	if(MCMC_PAR[i].Sampling_Switch==1)
		  	{
		  		if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
		  		{
		  			if(Sample_Physical_Parameters==1)
		  				if(Time_Dependant_PhysPar==1)
		  					for(snap=0;snap<NOUT;snap++)
		  						printf("p%d=%0.10g ",i+1,MCMC_PAR[i].PropValue[snap]);
		  				else
		  					printf("p%d=%0.10g ",i+1,MCMC_PAR[i].PropValue[0]);
		  			printf("\n");
		  		}
		  		else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
		  		{
		  			if(Sample_Cosmological_Parameters==1)
		  			{
		  				printf("cosm_p%d=%0.10g ",i+1,MCMC_PAR[i].PropValue[0]);
		  				printf("\n");
		  			}
		  		}
		  	}
		  printf("\n******************************************************\n\n\n\n");
		}
	}
	return (AcceptanceCounter);
}




///////////
//PROPOSE//
///////////


/*@brief Function to propose new parameters given by
 *       a random normal distribution gassdev*/

double propose_new_parameters ()
{
  double qratio;
  int i, snap;

  for(i=0;i<Npar;++i)
  	if(MCMC_PAR[i].Sampling_Switch==1)
  	{
  		if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
  		{
  			if(Sample_Physical_Parameters==0)
  				for(snap=0;snap<NOUT;snap++)
  					MCMC_PAR[i].PropValue[snap]= MCMC_PAR[i].Value[0];
  			else
  				for(snap=0;snap<NOUT;snap++)
  				{
  					if(Time_Dependant_PhysPar==0 && snap>0)
  						MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].PropValue[0];
  					else
  					{
  						do
  							MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap] * exp(MCMC_LogStep_Size * gassdev(&MCMCseed));
  						while(MCMC_PAR[i].PropValue[snap] < MCMC_PAR[i].PriorMin
  								|| MCMC_PAR[i].PropValue[snap] > MCMC_PAR[i].PriorMax);
  					}
  				}
  		}
  		else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
  			if(Sample_Cosmological_Parameters==0)
  				MCMC_PAR[i].PropValue[0]= MCMC_PAR[i].Value[0];
  	}

  if(Sample_Cosmological_Parameters==1)
  	do
  	{
  		 for(i=0;i<Npar;++i)
  			 if(MCMC_PAR[i].Sampling_Switch==1)
  				 if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
  					 do
  						 MCMC_PAR[i].PropValue[0] = MCMC_PAR[i].Value[0] * exp(MCMC_LogStep_Size * gassdev(&MCMCseed));
  					 while(MCMC_PAR[i].PropValue[0] < MCMC_PAR[i].PriorMin
  							 || MCMC_PAR[i].PropValue[0] > MCMC_PAR[i].PriorMax);

  		reset_cosmology();

  	}while(fabs(ZZ[ListOutputSnaps[0]]-ListOutputRedshifts[0]) > 0.1);


  //to make just some parameters the same at all z
  /*for(i=0;i<Npar;++i)
  {
  	for(snap=0;snap<NOUT;snap++)
  	{
  		//if(i<3 || (i>3 && i<6) || (i>6 && i<9) )
	    //if(i<6 || (i>6 && i<9) || i>9)
  		// if(i!=9)
  		prop[snap][i]=prop[0][i];
  	}
  }*/


  qratio = 1;
  //qratio= prop[0]/p[0]*prop[1]/p[1]*prop[2]/p[2]*prop[3]/p[3]*prop[4]/p[4];
  return qratio;
}

void initialize_mcmc_par ()
{
  int i, snap;
  double aux_p;
  char buf[1000];
  FILE *fa;

  //initialize structure to contain parameter names, values, priors and other properties
  MCMC_PAR = mymalloc("MCMC_PAR", sizeof(struct MCMC_PAR) * Npar);

  sprintf(buf, "%s", MCMCParameterValues);
  if(!(fa = fopen(buf, "r")))
  {
  	char sbuf[1000];
 	  sprintf(sbuf, "can't open file `%s'\n", buf);
 	  terminate(sbuf);
 	}

  fgets(buf, 300, fa);
  for(i=0;i<Npar;i++)
  {
  	fscanf(fa,"%s %lg %lg %lg %s %d\n",MCMC_PAR[i].Name, &MCMC_PAR[i].PropValue[0],
  				 &MCMC_PAR[i].PriorMin, &MCMC_PAR[i].PriorMax, MCMC_PAR[i].Type, &MCMC_PAR[i].Sampling_Switch);
  }
  fgets(buf, 300, fa);
  for(i=0;i<Npar;i++)
  	fscanf(fa,"%lg",&MCMC_PAR[i].Value[0]);

  fclose(fa);

  //if MCMC_Initial_Par_Displacement=1 introduce displacement in parameter values
  for(i=0;i<Npar;++i)
  	if(MCMC_PAR[i].Sampling_Switch==1)
  	{
  		if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
  			for(snap=0;snap<NOUT;snap++)
  			{
  				if(Sample_Physical_Parameters==1)
  				{
  					if(Time_Dependant_PhysPar==1 || snap==0)
  					{
  						aux_p=MCMC_PAR[i].Value[snap];
  						do
  							MCMC_PAR[i].Value[snap] = aux_p * exp(MCMC_Initial_Par_Displacement * gassdev(&MCMCseed));
  						while(MCMC_PAR[i].Value[snap] < MCMC_PAR[i].PriorMin
  								|| MCMC_PAR[i].Value[snap] > MCMC_PAR[i].PriorMax);
  						MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap];
  					}
  					else
  					{
  						MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];
  						MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[0];
  					}
  				}
  				else
  				{
  					MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];
  					MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[0];
  				}
  			}
  		else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
  			if(Sample_Cosmological_Parameters==1)
  			{
  				aux_p=MCMC_PAR[i].Value[snap];
  				do
  					MCMC_PAR[i].Value[snap] = aux_p * exp(MCMC_Initial_Par_Displacement * gassdev(&MCMCseed));
  				while(MCMC_PAR[i].Value[0] < MCMC_PAR[i].PriorMin
  						|| MCMC_PAR[i].Value[0] > MCMC_PAR[i].PriorMax);
  				MCMC_PAR[i].PropValue[0] = MCMC_PAR[i].Value[0];
  			}
  			else
  			{
  				MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];
  				MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[0];
  			}
}

  		//zlist_0012_0047.txt
   /* p[0][6]=6.6;
      p[1][6]=5.9;
      p[2][6]=3.7;
      p[3][6]=7.4;
      p[0][9]=0.99;
      p[1][9]=0.95;
      p[2][9]=0.6;
      p[3][9]=0.93;*/

}

void read_mcmc_par (int snapnum)
{
	int snap, i;

	//betwenn snapnum[i] and snapnum[i+1] parameters have the values of snap[i+1]
	for(snap=0;snap<NOUT;snap++)
		if(snapnum < ListOutputSnaps[NOUT-snap-1]+1)
			break;
	snap=NOUT-snap-1;

	if(snap<0)
		snap=0;

	for(i=0;i<Npar;i++)
		if(MCMC_PAR[i].Sampling_Switch==1)
		{
			if(strcmp(MCMC_PAR[i].Name,"SfrEfficiency")==0)
				SfrEfficiency = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"SfrBurstEfficiency")==0)
				SfrBurstEfficiency = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"SfrBurstSlope")==0)
				SfrBurstSlope = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"AgnEfficiency")==0)
			{
				AgnEfficiency = MCMC_PAR[i].PropValue[snap];
				AgnEfficiency /= (UnitMass_in_g / UnitTime_in_s * 1.58e-26);
			}
			else if(strcmp(MCMC_PAR[i].Name,"BlackHoleGrowthRate")==0)
				BlackHoleGrowthRate = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"BlackHoleCutoffVelocity")==0)
				BlackHoleCutoffVelocity = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"FeedbackReheatingEpsilon")==0)
				FeedbackReheatingEpsilon = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"ReheatPreVelocity")==0)
				ReheatPreVelocity = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"ReheatSlope")==0)
				ReheatSlope = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"FeedbackEjectionEfficiency")==0)
				FeedbackEjectionEfficiency = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"EjectPreVelocity")==0)
				EjectPreVelocity = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"EjectSlope")==0)
				EjectSlope = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"ReIncorporationFactor")==0)
				ReIncorporationFactor	= MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"ReincZpower")==0)
				ReincZpower = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"ReincVelocitypower")==0)
				ReincVelocitypower = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"Yield")==0)
				Yield = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"ThreshMajorMerger")==0)
				ThreshMajorMerger = MCMC_PAR[i].PropValue[snap];

			else if(strcmp(MCMC_PAR[i].Name,"Reionization_z0  ")==0)
				Reionization_z0 = MCMC_PAR[i].PropValue[snap];
			else if(strcmp(MCMC_PAR[i].Name,"Reionization_zr  ")==0)
				Reionization_zr = MCMC_PAR[i].PropValue[snap];
			//printf("EjectSlope=%g\n",EjectSlope);
	}
}



/*@brief Read in the IDs and Weights of the FOFs groups that
 *       constitute the sample for which galaxies will be
 *       compared with observations */
void read_sample_info (void)
{
  int DumbTreeNrColector, DumbFileNrColector;
  int i, snap;
  FILE *fa;
  char buf[1000];

  for(snap=0;snap<NOUT;snap++)
    {
      sprintf(buf, "%s/sample_allz_nh_%d%d.dat", MCMCTreeSampleDir, 
	      (MCMCTreeSampleFile-LastDarkMatterSnapShot)/100, ListOutputSnaps[snap]);
      if(!(fa = fopen(buf, "r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}
      
      fscanf(fa, "%d \n", &NTreesInSample[snap]);
      //only the sample for the last snap has all the trees
      if(NTreesInSample[snap]!=Ntrees && snap==LastDarkMatterSnapShot)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "NTreesInSample != Ntrees");
	  terminate(sbuf);
	}
      
      Weight[snap] = malloc(sizeof(double *) * (NTreesInSample[snap]));
      SampleIDs[snap] = malloc(sizeof(long long *) * (NTreesInSample[snap]));
      
      for(i=0;i<NTreesInSample[snap];i++)
	fscanf(fa, "%lld %d %d %lg\n", &SampleIDs[snap][i], &DumbTreeNrColector, &DumbFileNrColector, &Weight[snap][i]);
      
      fclose(fa);
    }
}



/*@brief Read in the arrays of observational data. They will be compared
 *       with the outputs from the SAM On the get_likelihood routine*/
void read_observations (void)
{
  int i, j, MaxBins, snap;
  FILE *f[NTests];
  double BinValueColector;
  char buf[1000];

  MaxBins = 0;
  //first just determine the observation with most bins, over allz, to allocate MCMC_Obs
  for(snap=0;snap<NOUT;snap++)
    {
	  //the round of ZZ[] is to assure that independently of the cosmology used you still
	  //round to z=1.0 or 2.0,etc...
	  for(i=0;i<NTests;++i)
	    {
		  sprintf(buf, "%s/observationaltest%d_z%1.2f.txt",ObsConstraintsDir,i,
				  (float)(int)(ZZ[ListOutputSnaps[snap]]+0.5) );
		  if((f[i]=fopen(buf,"r"))==NULL)
		    {
			  char sbuf[1000];
			  sprintf(sbuf, "can't open file `%s'\n", buf);
			  terminate(sbuf);
		    }
	    }
	  /* This values will give the size of the observational arrays
	   * They will be used in get_likelihood. MaxBins will be the
	   * largest value used to allocate memmory for the MCMC_Obs struct */
	  for(i = 0; i < NTests; ++i)
	    {
		  fscanf(f[i], "%d", &Nbins[snap][i]);
		  MaxBins = (Nbins[snap][i] > MaxBins ? Nbins[snap][i] : MaxBins);
		  fclose(f[i]);
	    }
    }

  //allocate structure to contain observational data (sizeof largest obs data set)
  MCMC_Obs = mymalloc("MCMC_Obs", sizeof(struct MCMC_OBSCONSTRAINTS) * MaxBins);

  //now read the observations
  for(snap=0;snap<NOUT;snap++)
     {
 	  //the round of ZZ[] is to assure that independently of the cosmology used you still
 	  //round to z=1.0 or 2.0,etc...
 	  for(i=0;i<NTests;++i)
 	    {
 		  sprintf(buf, "%s/observationaltest%d_z%1.2f.txt",ObsConstraintsDir,i,
 				  (float)(int)(ZZ[ListOutputSnaps[snap]]+0.5) );
 		  if((f[i]=fopen(buf,"r"))==NULL)
 		    {
 			  char sbuf[1000];
 			  sprintf(sbuf, "can't open file `%s'\n", buf);
 			  terminate(sbuf);
 		    }
 	    }

 	  /* This values will give the size of the observational arrays
 	   * They will be used in get_likelihood. MaxBins will be the
 	   * largest value used to allocate memmory for the MCMC_Obs struct */
 	  for(i = 0; i < NTests; ++i)
 		  fscanf(f[i], "%d", &Nbins[snap][i]);


	  for(i = 0; i < NTests; i++)
		OutputObs[snap][i] = malloc(sizeof(double *) * (ChainLength+1));

	  for(i = 0; i < NTests; i++)
		for(j = 0; j < ChainLength+1; j++)
		  OutputObs[snap][i][j] = malloc(sizeof(double *) * (MaxBins + 1));


	  //Read observational data
	  for(i = 0; i < NTests; i++)
	    {
		  if(i < NChiTests)//Chi2 Tests
	    	for(j = 0; j < Nbins[snap][i]; j++)
	   		  fscanf(f[i], "%lg %lg %lg", &MCMC_Obs[j].Bin[snap][i], &MCMC_Obs[j].Obs[snap][i], &MCMC_Obs[j].Error[snap][i]);

		  else  //Binomial tests
			for(j = 0; j < Nbins[snap][i]; j++)
		      fscanf(f[i], "%lg %lg %lg", &BinValueColector, &MCMC_Obs[j].ObsUp[snap][i - NChiTests],
		    		  &MCMC_Obs[j].ObsDown[snap][i - NChiTests]);

		  fclose(f[i]);
	    }
    }
}



/*For each step (for each parameter step), print the semi-analytic predictions
  used for the tests*/
void write_comparison_to_obs()
{
	int i, j, k, snap;
	FILE *fa;
	char buf[1000];

	sprintf(buf, "%sobs_output_%d.txt", OutputDir, ThisTask);
	if((fa = fopen(buf, "w")) == NULL)
	{
		char sbuf[1000];
		sprintf(sbuf, "can't open file `%s'\n", buf);
		terminate(sbuf);
    }

	fprintf(fa, " %d %d\n", NTests, ChainLength+1);
	for(i = 0; i < NTests; i++)
	  for(snap=0;snap<NOUT;snap++)
		fprintf(fa, " %d\n", Nbins[snap][i]);

	for(snap=0;snap<NOUT;snap++)
	  for(i = 0; i < NTests; i++)
		for(j = 0; j < ChainLength+1; j++)
		  {
			for(k = 0; k < Nbins[snap][i]; k++)
				fprintf(fa, " %g", OutputObs[snap][i][j][k]);
			fprintf(fa, " %g\n", OutputObs[snap][i][j][k]);
		  }

	for(snap=0;snap<NOUT;snap++)
	  for(i = 0; i < NTests; i++)
		free(OutputObs[snap][i]);


	fclose(fa);
}




//////////
//GASDEV//
//////////

//Gives a random normal deviate using ran3 (ran1 NR)


double gassdev(long *idum)
{
  static int iset = 0;
  static double gset;
  double fac, r, v1, v2;

  if(iset == 0)
    {
      do
	{
	  v1 = 2.0 * ran3(idum) - 1.0;
	  v2 = 2.0 * ran3(idum) - 1.0;
	  r = v1 * v1 + v2 * v2;
	}
      while(r >= 1.0 || r == 0.0);
      fac = sqrt(-2.0 * log(r) / r);

      //Box Muller deviates to get two normal deviates
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }

  else
    {
      iset = 0;
      return gset;
    }
}






////////
//RAN3//
////////

double ran3(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
	*idum = 1;
      else
	*idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ;
	  *idum = IA * (*idum - k * IQ) - IR * k;
	  if(*idum < 0)
	    *idum += IM;
	  if(j < NTAB)
	    iv[j] = *idum;
	}

      iy = iv[0];
    }

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;

  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}
#endif //MCMC



#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef RNMX
#undef EPS
