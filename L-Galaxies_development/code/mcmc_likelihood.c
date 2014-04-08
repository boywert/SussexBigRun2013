#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"
#include "mcmc_halomodel.h"

#define FPMIN 1.0e-30
#define ASWITCH 100
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAXIT 100




//////////////
//LIKELIHOOD//
//////////////


/** @file mcmc_likelihood.c
 *  @brief This function computes the likelihhod of SAM galaxies in respect
 *         to a number of observational constraints (Henriques2009)
 *
 *  This function computes the likelihood between the SAM and a set of observational
 *  data sets. The galaxy properties for the SAM are computed using a given set
 *  of input parameters and in the end of each run this function is called to get
 *  the likelihood of properties with the new set of parameters and compare it with
 *  the previous run.
 *
 * Three different tests are available: Chi^2, Maximum Likelihhod Method and Binomial
 * probability. These are used to compare the model with the different observational
 * properties read at read_observations and stored in struct MCMC_Obs.
 *
 * prob[0:13] are computed and the likelihood is given by a desired
 * combination of these at the exit of the get_likelihho() function.
 *
 * 1st Test - Chisquare test for the SMF
 * MCMC_Obs[].Obs[0] & MCMC_Obs[].Error[0] - observations from baldry2008
 * MCMC_GAL[].StellarMass - Masses from SAM
 *
 * 8th Test - MLM test for the colours
 * MCMC_Obs[].Obs[1] & MCMC_Obs[].Error[1] - observations from baldry2004
 * binredfraction - fraction of red galaxies in each mass bin
 * fracterr - error in the SAM red fraction (assumed to have the same value in all bins
 * - 0.05 to account for the wiggles in the 2 color luminosity function, due to the fact
 * that only one file is being used )
 *
 * 2nd Test - Chisquare test for the K-Band LF
 * MCMC_Obs[].Obs[2] & MCMC_Obs[].Error[2] - observations from Cole2003 + Bell2003+ Jones2006
 * MCMC_GAL[].MagK - SAM MagK
 *
 *
 * 10th Test - Binomial test for the Black Hole-Bulge Mass relation
 * MCMC_Obs[].ObsUp[0] & MCMC_Obs[].ObsDown[0] - observations for the BHBM from Haring & Rix 2004
 * MCMC_GAL[].Bulge & MCMC_GAL[].BlackHoleMass - masses from SAM
 * binblackholeup - numbers of galaxies in the two upper bins on the bulge-blackhole mass relation
 * binblackholedown - numbers of galaxies in the two lower bins on the bulge-blackhole mass relation */

#ifdef MCMC
double get_likelihood()
{
	//Variables for all the LF/SM function like tests
	double *binsamdata, *samdata;
	//variables for the bulge-blackhole mass test using binomial
	double binblackholeup[2]={0.0,0.0}, binblackholedown[2]={0.0,0.0};
	double final_probability, redshift_probability, current_probability;
	double prob_SMF_z0;
	int i, j, snap;
	double color_UV, color_VJ;
	//OBSERVATIONAL CUT
	//double offset_color_cut[NOUT]={0.00, 0.69, 0.59, 0.59, 0.59}; //not used at z=0
	//double slope_color_cut[NOUT]={0.00, 0.88, 0.88, 0.88, 0.88};
	//BEST FIT TO MODEL CUT
#ifndef HALOMODEL
	double offset_color_cut[NOUT]={0.00, 1.085, 1.1, 1.0, 1.15}; //not used at z=0
	double slope_color_cut[NOUT]={0.00, 0.5, 0.48, 0.38, 0.18};
#else
	double offset_color_cut[NOUT]={0.00}; //not used at z=0
	double slope_color_cut[NOUT]={0.00};
#endif

	/* Bin samdata into binsamdata according to observational constraints.
	 * The first argument of the bin functions and of the likelihood
	 * function (Chi^2 ot MLM) indicates the observational data set to use.*/

	//ALL Luminosity Function to be compared in units of M-5logh and phi(Mpc-3h-3mag-1)

	final_probability=1.;
	for(snap=0;snap<NOUT;snap++)
	{
		printf("snap=%d\n",snap);
		redshift_probability=1.;

		if((samdata = malloc(sizeof(double) * TotMCMCGals[snap])) == NULL)
			terminate("get_likelihood");

#ifdef HALOMODEL
		correct_for_correlation(snap);
#endif

		for(i=0;i<MCMCNConstraints;i++)
		{
			if(MCMC_Obs[i].ObsTest_Switch_z[snap]==1)
			{

				if((binsamdata = malloc(sizeof(double) * Nbins[snap][i])) == NULL)
					terminate("get_likelihood");


/******************************
 **     Chi_Sq TESTS         **
 ******************************/
				if(strcmp(MCMC_Obs[i].TestType,"chi_sq")==0)
				{

					//bin all the galaxie into binsamdata for properties that just require normal histograms
					for(j = 0; j < TotMCMCGals[snap]; j++)
					{
						samdata[j] = 0.; //initialize
						if(strcmp(MCMC_Obs[i].Name,"StellarMassFunction")==0)
							samdata[j] = MCMC_GAL[j].StellarMass[snap];
						if(strcmp(MCMC_Obs[i].Name,"KBandLF")==0)
							samdata[j] = MCMC_GAL[j].MagK[snap];
						if(strcmp(MCMC_Obs[i].Name,"BBandLF")==0)
						{
							if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)<0.2) //Bj band at z=0
								samdata[j]=MCMC_GAL[j].MagB[snap]-0.28*(MCMC_GAL[j].MagB[snap]-MCMC_GAL[j].MagV[snap]);
							else //BBand at all other z
								samdata[j]=MCMC_GAL[j].MagB[snap];
						}
						if(strcmp(MCMC_Obs[i].Name,"uBandLF")==0)
							samdata[j] = MCMC_GAL[i].Magu[snap];
						if(strcmp(MCMC_Obs[i].Name,"ColdGasMassFunction")==0)
							samdata[j] = MCMC_GAL[j].ColdGas[snap];

						//Stellar Mass Function of Passive Galaxies
						if(strcmp(MCMC_Obs[i].Name,"StellarMassFunctionPassive")==0)
								if( MCMC_GAL[j].Sfr[snap]* 1.e9 /pow(10.,MCMC_GAL[j].StellarMass[snap]) < 0.01 )
									samdata[j] = MCMC_GAL[j].StellarMass[snap];

						//Stellar Mass Function of Active Galaxies
						if(strcmp(MCMC_Obs[i].Name,"StellarMassFunctionActive")==0)
								if( MCMC_GAL[j].Sfr[snap]* 1.e9 /pow(10.,MCMC_GAL[j].StellarMass[snap]) > 0.3 )
									samdata[j] = MCMC_GAL[j].StellarMass[snap];

						//Stellar Mass Function of Red Galaxies
						//original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[j].Magr[snap]+20.07)/1.09))
						if(strcmp(MCMC_Obs[i].Name,"StellarMassFunctionRed")==0)
						{
							if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)<0.2) //z=0
							{
								if( (MCMC_GAL[j].Magu[snap]-MCMC_GAL[j].Magr[snap]) > (1.9-0.244*tanh((MCMC_GAL[j].Magr[snap]+20.07)/1.09)))
									samdata[j] = MCMC_GAL[j].StellarMass[snap];
							}
							else //z>0
							{
								color_UV=(MCMC_GAL[j].MagU[snap]-MCMC_GAL[j].MagV[snap]);
								color_VJ=(MCMC_GAL[j].MagV[snap]-MCMC_GAL[j].MagJ[snap]);
								if( (color_VJ < (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV > 1.3) ||
										(color_VJ > (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV > color_VJ*slope_color_cut[snap]+offset_color_cut[snap]) )
									samdata[j] = MCMC_GAL[j].StellarMass[snap];
							}
						}

						//Stellar Mass Function of Blue Galaxies
						//original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[j].Magr[snap]+20.07)/1.09))
						if(strcmp(MCMC_Obs[i].Name,"StellarMassFunctionBlue")==0)
						{
							if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)<0.2) //z=0
							{
								if( (MCMC_GAL[j].Magu[snap]-MCMC_GAL[j].Magr[snap]) < (1.9-0.244*tanh((MCMC_GAL[j].Magr[snap]+20.07)/1.09)))
									samdata[j] = MCMC_GAL[j].StellarMass[snap];
							}
							else //z>0
							{
								color_UV=(MCMC_GAL[j].MagU[snap]-MCMC_GAL[j].MagV[snap]);
								color_VJ=(MCMC_GAL[j].MagV[snap]-MCMC_GAL[j].MagJ[snap]);

								if( (color_VJ < (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV < 1.3) ||
										(color_VJ > (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV < color_VJ*slope_color_cut[snap]+offset_color_cut[snap]) )
									samdata[j] = MCMC_GAL[j].StellarMass[snap];
							}
						}

					}//end loop on number of galaxies
					bin_function(i, binsamdata, samdata, snap);

					//SFRD only has one bin
					if(strcmp(MCMC_Obs[i].Name,"SFRD")==0)
					{
						binsamdata[0]=0;
						for(j = 0; j < TotMCMCGals[snap]; j++)
							binsamdata[0]+=MCMC_GAL[j].Sfr[snap]*MCMC_GAL[j].Weight[snap];
					}

#ifdef HALOMODEL
					//Correlation Function - requires more than just binning
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_8.77_9.27",28)==0)
						compute_correlation_func(i, binsamdata, snap, 8.77, 9.27);
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_9.27_9.77",28)==0)
						compute_correlation_func(i, binsamdata, snap, 9.27, 9.77);
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_9.77_10.27",28)==0)
						compute_correlation_func(i, binsamdata, snap, 9.77, 10.27);
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_10.27_10.77",28)==0)
						compute_correlation_func(i, binsamdata, snap, 10.27, 10.77);
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_10.77_11.27",28)==0)
						compute_correlation_func(i, binsamdata, snap, 10.77, 11.27);
					if(strncmp(MCMC_Obs[i].Name,"Clustering_MassBins_11.27_11.47",28)==0)
						compute_correlation_func(i, binsamdata, snap, 11.27, 11.77);
#endif

					current_probability=chi_square_probability(i, binsamdata, snap);

				}//end chi_sq tests


/******************************
 ** Maximum Likelihood TESTS **
 ******************************/
				if(strcmp(MCMC_Obs[i].TestType,"maxlike")==0)
				{
					if(strcmp(MCMC_Obs[i].Name,"RedFraction")==0)
					{
						bin_red_fraction(i, binsamdata, snap);
						current_probability=maximum_likelihood_probability(i, binsamdata, snap);
					}
				}

				if(strcmp(MCMC_Obs[i].TestType,"maxlike")==0)
				{
					if(strcmp(MCMC_Obs[i].Name,"PassiveFraction")==0)
					{
						bin_passive_fraction(i, binsamdata, snap);
						current_probability=maximum_likelihood_probability(i, binsamdata, snap);
					}
				}

/******************************
 **    Binomial TESTS        **
 ******************************/
				if(strcmp(MCMC_Obs[i].TestType,"binomial")==0)
				{
					if(strcmp(MCMC_Obs[i].Name,"BlackHoleBulgeMass")==0)
					{
						bin_bhbm(binblackholeup, binblackholedown, snap);
						current_probability=binomial_probability(i, binblackholeup, binblackholedown, snap);
					}
				}
//END BINOMIAL TESTS


				free(binsamdata);


				redshift_probability*=current_probability;
				printf("prob[%d]=%0.5e\n",i,current_probability);
				//write likelihood for each constraint at each MCMC step
				fprintf(FILE_MCMC_PredictionsPerStep[snap][i], " %0.5e", current_probability);

				//save probability of the stellar mass function at z=0
				//if(snap==0 && strcmp(MCMC_Obs[i].Name,"StellarMassFunction")==0)
				//	prob_SMF_z0=current_probability;

			}//end if(MCMC_Obs[i].ObsTest_Switch_z[snap]==1)
		}//end loop on MCMCNConstraints

		printf("snap[%d] prob=%0.5e\n",snap, redshift_probability);

		final_probability*=redshift_probability;
		free(samdata);
#ifdef HALOMODEL
		free(MCMC_FOF2);
		free(HashTable);
#endif
	}//end loop on snaps


	printf("final prob=%0.5e\nlog_like=%0.8g\n",final_probability, -log10(final_probability));

	//write total likelihood into separate files and into
	//the end of the line on the file with the comparison to each constraint
	fprintf(FILE_MCMC_LIKELIHOOD,"%0.8g\n",-log10(final_probability));
	fflush(FILE_MCMC_LIKELIHOOD);
	for(i=0;i<MCMCNConstraints;i++)
		for(snap=0;snap<NOUT;snap++)
			if(MCMC_Obs[i].ObsTest_Switch_z[snap]==1)
			{
				fprintf(FILE_MCMC_PredictionsPerStep[snap][i], " %0.8g\n", -log10(final_probability));
				fflush(FILE_MCMC_PredictionsPerStep[snap][i]);
			}
	//if(final_probability>1.e-22 && prob_SMF_z0>1.e-4)
	//	printf("ola\n");

	return final_probability;
}

#ifdef HALOMODEL
void correct_for_correlation(int snap)
{
	int jj, *CumulativeNgals;
	int fof, currentNfof=0, TotalNFof=0;

	if((HashTable = malloc(sizeof(int) * TotMCMCGals[snap])) == NULL)
		terminate("correct_for_correlation");

//	printf("Nfofs=%d\n",NFofsInSample[snap]);
	for(fof=0;fof<NFofsInSample[snap]; fof++)
		if(MCMC_FOF[fof].NGalsInFoF[snap]>0)
			TotalNFof++;

	MCMC_FOF2 = malloc(sizeof(struct MCMC_FOF_struct) * TotalNFof);
	CumulativeNgals = malloc(sizeof(int) * TotalNFof);

	currentNfof=0;
	CumulativeNgals[0]=0;

	//Assign FOF properties to a different structure, in case not all FOFs have been computed
	//If not all the trees were selected in main.c
	for(fof=0;fof<NFofsInSample[snap]; fof++)
		if(MCMC_FOF[fof].NGalsInFoF[snap]>0)
		{
			if(currentNfof<TotalNFof-1)
				CumulativeNgals[currentNfof+1]=CumulativeNgals[currentNfof]+MCMC_FOF[fof].NGalsInFoF[snap];

			MCMC_FOF2[currentNfof].NGalsInFoF[snap]=	MCMC_FOF[fof].NGalsInFoF[snap];
			MCMC_FOF2[currentNfof].M_Crit200[snap]=	MCMC_FOF[fof].M_Crit200[snap];
			MCMC_FOF2[currentNfof].M_Mean200[snap]=	MCMC_FOF[fof].M_Mean200[snap];
			MCMC_FOF2[currentNfof].IndexOfCentralGal[snap]=	MCMC_FOF[fof].IndexOfCentralGal[snap];
			//printf("FOF[%d], Ngals =%d\n",currentNfof,MCMC_FOF2[currentNfof].NGalsInFoF[snap]);
			for(jj=0;jj<TotMCMCGals[snap];jj++)
				if(MCMC_GAL[jj].fofid[snap] == fof)
				{
					MCMC_GAL[jj].fofid[snap] = currentNfof;
					MCMC_GAL[jj].ngal[snap] = MCMC_FOF[fof].NGalsInFoF[snap];
					MCMC_GAL[jj].M_Crit200[snap] = MCMC_FOF[fof].M_Crit200[snap];
					MCMC_GAL[jj].M_Mean200[snap] = MCMC_FOF[fof].M_Mean200[snap];
				}
			currentNfof++;
		}

	UsedFofsInSample[snap]=TotalNFof;

	//BUILD HASHTABLE

	for(jj=0;jj<TotMCMCGals[snap];jj++)
		HashTable[jj]=-1;

	int kk;

	for(fof=0;fof<currentNfof; fof++)
	//	for(fof=0;fof<1; fof++)
	{
			//printf("CumGals=%d\n",CumulativeNgals[fof]);
		for(jj=0;jj<TotMCMCGals[snap];jj++)
		{
			if(MCMC_GAL[jj].fofid[snap]==fof)
			{
			  //if type=0 it gets the first place in hashtable for this group (CumulativeNgals[fof])
				if(MCMC_GAL[jj].Type[snap]==0)
				{
					HashTable[CumulativeNgals[fof]]=jj;
					MCMC_FOF2[fof].IndexOfCentralGal[snap]=CumulativeNgals[fof];
					//printf("fof=%d indexofcentral=%d\n",fof, MCMC_FOF2[fof].IndexOfCentralGal[snap]);
				}
			//if not gets the first available place in hashtable for this group (CumulativeNgals[fof]+kk)
				else
					for(kk=1;kk<MCMC_FOF2[fof].NGalsInFoF[snap];kk++)
						if(HashTable[CumulativeNgals[fof]+kk]==-1)
						{
							HashTable[CumulativeNgals[fof]+kk]=jj;
							break;
						}
			//printf("Gal[%d] FofId =%d Gal_Type=%d Ngals in group=%d\n",jj,MCMC_GAL[jj].fofid[snap],MCMC_GAL[jj].Type[snap],MCMC_GAL[jj].ngal[snap]);
			}
		}
	}//loop on fof to get hashtable


	free(CumulativeNgals);
/*
	FILE *fa;
	char buf[1000];

	//Write MCMC gal structure
	sprintf(buf, "%s/MCMC_GAL_z0", OutputDir);
	if((fa = fopen(buf, "w")) == NULL)
	{
		char sbuf[1000];
		sprintf(sbuf, "can't open file `%s'\n", buf);
		terminate(sbuf);
	}

	//write
	MCMC_GAL2 = mymalloc("MCMC_Gal_reorder", sizeof(struct MCMC_GALAXY_reorder) * TotMCMCGals[snap]);

	for(jj=0;jj<TotMCMCGals[snap];jj++)
	{
		MCMC_GAL2[jj].StellarMass[snap]=MCMC_GAL[HashTable[jj]].StellarMass[snap];
		MCMC_GAL2[jj].ColdGas[snap]=MCMC_GAL[HashTable[jj]].ColdGas[snap];
		MCMC_GAL2[jj].BulgeMass[snap]=MCMC_GAL[HashTable[jj]].BulgeMass[snap];
		MCMC_GAL2[jj].BlackHoleMass[snap]=MCMC_GAL[HashTable[jj]].BlackHoleMass[snap];
		MCMC_GAL2[jj].Sfr[snap]=MCMC_GAL[HashTable[jj]].Sfr[snap];
		MCMC_GAL2[jj].MagU[snap]=MCMC_GAL[HashTable[jj]].MagU[snap];
		MCMC_GAL2[jj].MagB[snap]=MCMC_GAL[HashTable[jj]].MagB[snap];
		MCMC_GAL2[jj].MagV[snap]=MCMC_GAL[HashTable[jj]].MagV[snap];
		MCMC_GAL2[jj].MagJ[snap]=MCMC_GAL[HashTable[jj]].MagJ[snap];
		MCMC_GAL2[jj].MagK[snap]=MCMC_GAL[HashTable[jj]].MagK[snap];
		MCMC_GAL2[jj].Magu[snap]=MCMC_GAL[HashTable[jj]].Magu[snap];
		MCMC_GAL2[jj].Magg[snap]=MCMC_GAL[HashTable[jj]].Magg[snap];
		MCMC_GAL2[jj].Magr[snap]=MCMC_GAL[HashTable[jj]].Magr[snap];
		MCMC_GAL2[jj].Magi[snap]=MCMC_GAL[HashTable[jj]].Magi[snap];
		MCMC_GAL2[jj].Magz[snap]=MCMC_GAL[HashTable[jj]].Magz[snap];
		MCMC_GAL2[jj].Weight[snap]=MCMC_GAL[HashTable[jj]].Weight[snap];
		MCMC_GAL2[jj].fofid[snap]=MCMC_GAL[HashTable[jj]].fofid[snap];
		MCMC_GAL2[jj].M_Crit200[snap]=MCMC_GAL[HashTable[jj]].M_Crit200[snap];
		MCMC_GAL2[jj].x[snap]=MCMC_GAL[HashTable[jj]].x[snap];
		MCMC_GAL2[jj].y[snap]=MCMC_GAL[HashTable[jj]].y[snap];
		MCMC_GAL2[jj].z[snap]=MCMC_GAL[HashTable[jj]].z[snap];
		MCMC_GAL2[jj].Type[snap]=MCMC_GAL[HashTable[jj]].Type[snap];
		MCMC_GAL2[jj].ngal[snap]=MCMC_GAL[HashTable[jj]].ngal[snap];
	}

	fwrite(&TotMCMCGals[snap], sizeof(int), 1, fa);
	fwrite(MCMC_GAL2, TotMCMCGals[snap], sizeof(struct MCMC_GALAXY_reorder), fa);

	myfree(MCMC_GAL2);

	fclose(fa);

	//Write FOF structure
	sprintf(buf, "%s/MCMC_FOF_z0", OutputDir);
	if((fa = fopen(buf, "w")) == NULL)
	{
		char sbuf[1000];
		sprintf(sbuf, "can't open file `%s'\n", buf);
		terminate(sbuf);
	}

	fwrite(&NFofsInSample[snap], sizeof(int), 1, fa);
	fwrite(MCMC_FOF2, NFofsInSample[snap], sizeof(struct MCMC_FOF_struct), fa);

	fclose(fa);
*/
}

void compute_correlation_func(int ObsNr, double *binsamdata, int snap, float mingalmass, float maxgalmass)
{
	double *r,*proj,*r_tmp,*proj_tmp;
	int ibin, ii, jj;
	char buf[1000];
	FILE *fa;
	gsl_spline *Proj_Spline;
	gsl_interp_accel *Proj_SplineAcc;

	NR=60;

	//printf("\ncalculating correlation function %0.2f < M < %0.2f\n",mingalmass,maxgalmass);

	r=malloc(NR*sizeof(double));
	proj=malloc(NR*sizeof(double));

	halomodel(r,proj,mingalmass,maxgalmass,snap);

	Proj_Spline=gsl_spline_alloc(gsl_interp_cspline,NR);
	Proj_SplineAcc=gsl_interp_accel_alloc();
	gsl_spline_init(Proj_Spline,r,proj,NR);

	//for(ii=0;ii<Nbins[snap][ObsNr]-1;ii++)
	//	printf("Nbins=%d r=%g\n",Nbins[snap][ObsNr], MCMC_Obs[ObsNr].Bin_low[snap][ii]+(MCMC_Obs[ObsNr].Bin_high[snap][ii]-MCMC_Obs[ObsNr].Bin_low[snap][ii])/2.)

	 for(ii=0;ii<Nbins[snap][ObsNr]-1;ii++)
		 binsamdata[ii]=gsl_spline_eval(Proj_Spline,MCMC_Obs[ObsNr].Bin_low[snap][ii]+(MCMC_Obs[ObsNr].Bin_high[snap][ii]-MCMC_Obs[ObsNr].Bin_low[snap][ii])/2.,Proj_SplineAcc);

#ifndef PARALLEL
	 //full3 - full PLANCK
	sprintf(buf, "%s/correlation_guo10_nfull4_%0.2f_%0.2f.txt",OutputDir, mingalmass,maxgalmass);
	if(!(fa = fopen(buf, "w")))
	{
		char sbuf[1000];
		sprintf(sbuf, "can't open file `%s'\n", buf);
		terminate(sbuf);
	}
	for(ii=0;ii<Nbins[snap][ObsNr]-1;ii++)
		fprintf(fa, "%g %g %g\n", MCMC_Obs[ObsNr].Bin_low[snap][ii]+(MCMC_Obs[ObsNr].Bin_high[snap][ii]-MCMC_Obs[ObsNr].Bin_low[snap][ii])/2.,binsamdata[ii],binsamdata[ii]*0.1);
	//for(ii=0;ii<NR;ii++)
	//fprintf(fa, "%g %g %g\n", r[ii],proj[ii],proj[ii]*0.1);
	fclose(fa);
#endif

	//original wrp calculation out of halo_model
	//printf("\n\n %g<M<%g",mingalmass,maxgalmass);
	//for(ii=0;ii<NR;ii++)
	//	printf("r=%g proj=%g\n",r[ii],proj[ii]);
  //interpolated into obs bins
	//for(ii=0;ii<Nbins[snap][ObsNr]-1;ii++)
	//	printf("%g %g %g\n", MCMC_Obs[ObsNr].Bin_low[snap][ii]+(MCMC_Obs[ObsNr].Bin_high[snap][ii]-MCMC_Obs[ObsNr].Bin_low[snap][ii])/2.,binsamdata[ii],binsamdata[ii]*0.1);

	 free(r);
	 free(proj);
	 gsl_spline_free(Proj_Spline);
	 gsl_interp_accel_free(Proj_SplineAcc);
}
#endif




/**@brief Bin Luminosity or Stellar Mass Function*/
void bin_function(int ObsNr, double *binsamdata, double *samdata, int snap)
{
  int i, k;
  FILE *fa;
  char buf[1000];

#ifndef PARALLEL
  sprintf(buf, "%s/mcmc_plus_obs%d_z%1.2f.txt",OutputDir,ObsNr,(double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.));
 	  if(!(fa = fopen(buf, "w")))
 	    {
 		  char sbuf[1000];
 		  sprintf(sbuf, "can't open file `%s'\n", buf);
 		  terminate(sbuf);
 	    }
#endif

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      binsamdata[i] = 0;

      for(k=0;k<TotMCMCGals[snap];k++)
    	  	if(samdata[k]>=MCMC_Obs[ObsNr].Bin_low[snap][i] && samdata[k] <= MCMC_Obs[ObsNr].Bin_high[snap][i])
    	  		binsamdata[i]+=MCMC_GAL[k].Weight[snap];
      //binsamdata[i]/=(BoxSize*BoxSize*BoxSize*Bin);
      binsamdata[i]/=(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i]);

    	fprintf(FILE_MCMC_PredictionsPerStep[snap][ObsNr], " %0.5e", binsamdata[i]);

#ifndef PARALLEL
      fprintf(fa, "%g %g %g %g\n", MCMC_Obs[ObsNr].Bin_low[snap][i]+(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i],
    		  	                   MCMC_Obs[ObsNr].Error[snap][i], binsamdata[i]);
#endif

      //printf("%g %g %g %g\n", MCMC_Obs[ObsNr].Bin_low[snap][i]+(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i],
      //	  	                   MCMC_Obs[ObsNr].Error[snap][i], binsamdata[i]);
    }
#ifndef PARALLEL
  fclose(fa);
#endif
}

/**@brief Bin fraction of red galaxies.
 *    The separation is done based on g-r color
 *    at z=0 and and U-V vs V-J at higher z */
void bin_red_fraction(int ObsNr, double *binredfraction, int snap)
{
  int i, k, IsRedGalaxy;
  double red, blue;
  double color, color_UV, color_VJ;
  //OBSERVATIONAL CUT
  //double offset_color_cut[snap]={0.00, 0.69, 0.59, 0.59, 0.59}; //not used at z=0
  //double slope_color_cut[snap]={0.00, 0.88, 0.88, 0.88, 0.88};
  //BEST FIT TO MODEL CUT
#ifndef HALOMODEL
	double offset_color_cut[NOUT]={0.00, 1.085, 1.1, 1.0, 1.15}; //not used at z=0
	double slope_color_cut[NOUT]={0.00, 0.5, 0.48, 0.38, 0.18};
#else
	double offset_color_cut[NOUT]={0.00}; //not used at z=0
	double slope_color_cut[NOUT]={0.00};
#endif

  FILE *fa;
  char buf[1000];

#ifndef PARALLEL
  sprintf(buf, "%s/mcmc_plus_obs%d_z%1.2f.txt",OutputDir,ObsNr,(double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.));
 	  if(!(fa = fopen(buf, "w")))
 	    {
 		  char sbuf[1000];
 		  sprintf(sbuf, "can't open file `%s'\n", buf);
 		  terminate(sbuf);
 	    }
#endif

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      binredfraction[i] = 0.;
      red = 0.;
      blue = 0.;

      for(k=0;k<TotMCMCGals[snap];k++)
      {
      	//color g-r cut for z=0  - baldry 2004
      	//original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[j].Magr[snap]+20.07)/1.09))
      	if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)<0.2)
      	{
      		color=MCMC_GAL[k].Magu[snap]-MCMC_GAL[k].Magr[snap];
      		IsRedGalaxy=(color > (1.9-0.244*tanh((MCMC_GAL[k].Magr[snap]+20.07)/1.09)));
      	}
      	//U-V vs V-J at higher redshift
      	else
      	{
      		color_UV=(MCMC_GAL[k].MagU[snap]-MCMC_GAL[k].MagV[snap]);
      		color_VJ=(MCMC_GAL[k].MagV[snap]-MCMC_GAL[k].MagJ[snap]);
      		IsRedGalaxy=( (color_VJ < (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV > 1.3 ) ||
      				(color_VJ > (1.3-offset_color_cut[snap])/slope_color_cut[snap] && color_UV > color_VJ*slope_color_cut[snap]+offset_color_cut[snap]) );
      	}

      	if(IsRedGalaxy == 1)
      	{
      		if(MCMC_GAL[k].StellarMass[snap]>=MCMC_Obs[ObsNr].Bin_low[snap][i]
      		   && MCMC_GAL[k].StellarMass[snap] <= MCMC_Obs[ObsNr].Bin_high[snap][i]) red+=MCMC_GAL[k].Weight[snap];
      	}
      	else
      	{
      		if(MCMC_GAL[k].StellarMass[snap]>=MCMC_Obs[ObsNr].Bin_low[snap][i]
      				&& MCMC_GAL[k].StellarMass[snap] <= MCMC_Obs[ObsNr].Bin_high[snap][i]) blue+=MCMC_GAL[k].Weight[snap];
      	}
      }

      if((blue + red) > 0)
      	binredfraction[i] = red * 1.0 / (blue * 1.0 + red * 1.0);
      else
      	binredfraction[i] = 0.;
      //printf("mass=%f obs =%f samcolors=%f\n", (MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i], binredfraction[i]);

    	fprintf(FILE_MCMC_PredictionsPerStep[snap][ObsNr], " %0.5e", binredfraction[i]);


#ifndef PARALLEL
      fprintf(fa, "%g %g %g %g\n",  MCMC_Obs[ObsNr].Bin_low[snap][i]+(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i],
    		  	                   MCMC_Obs[ObsNr].Error[snap][i], binredfraction[i]);
#endif

    }

#ifndef PARALLEL
  fclose(fa);
#endif

}

/**@brief Bin fraction of passive galaxies.
 *    The separation is done based on sfr*/
void bin_passive_fraction(int ObsNr, double *binpassivefraction, int snap)
{
  int i, k, IsPassiveGalaxy;
  double passive, active;
  double color;
  //if (strcmp(cq,"Cold")==0) {

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      binpassivefraction[i] = 0.;
      passive = 0.;
      active = 0.;

      for(k=0;k<TotMCMCGals[snap];k++)
      {
      	IsPassiveGalaxy=(MCMC_GAL[k].Sfr[snap]* 1.e9 /pow(10.,MCMC_GAL[k].StellarMass[snap]) < 0.01);

      	if(IsPassiveGalaxy == 1)
    	  		{
    	    		if(MCMC_GAL[k].StellarMass[snap]>=MCMC_Obs[ObsNr].Bin_low[snap][i]
    	    			&& MCMC_GAL[k].StellarMass[snap] <= MCMC_Obs[ObsNr].Bin_high[snap][i]) passive += MCMC_GAL[k].Weight[snap];
    	  		}
    	    else
    	    	{
    	    		if(MCMC_GAL[k].StellarMass[snap]>=MCMC_Obs[ObsNr].Bin_low[snap][i]
    	    			&& MCMC_GAL[k].StellarMass[snap] <= MCMC_Obs[ObsNr].Bin_high[snap][i]) active += MCMC_GAL[k].Weight[snap];
    	    	}
      }

      if((passive + active) > 0)
      	binpassivefraction[i] = passive * 1.0 / (active * 1.0 + passive * 1.0);
      else
      	binpassivefraction[i] = 0.;
      //printf("mass=%f obs =%f samcolors=%f\n",(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i], binpassivefraction[i]);
    }
}

/**@brief Bin SAM colours according to observations of Baldry2004*/
void bin_color_hist(int ObsNr, double *bincolorhist, int snap)
{
  int i, k;
  double color;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      bincolorhist[i] = 0;

      for(k=0;k<TotMCMCGals[snap];k++)
      {
    	  color=MCMC_GAL[k].Magg[snap]-MCMC_GAL[k].Magr[snap];
    	  if(color>=MCMC_Obs[ObsNr].Bin_low[snap][i]
    			  && color <= MCMC_Obs[ObsNr].Bin_high[snap][i]
    			  && MCMC_GAL[k].StellarMass[snap] > 9.0 && MCMC_GAL[k].StellarMass[snap] < 9.5)
    		  bincolorhist[i]=bincolorhist[i]+MCMC_GAL[k].Weight[snap];
      }
      //printf("%f %f %f\n",(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2., MCMC_Obs[ObsNr].Obs[snap][i],bincolorhist[i]);
    }

}

/**@brief Bin SAM BH and BM according to observations of Haring and Rix 2004*/
void bin_bhbm(double *binblackholeup, double *binblackholedown, int snap)
{
  int k;

  for(k = 0; k < TotMCMCGals[snap]; k++)
    {
      if(MCMC_GAL[k].BlackHoleMass[snap]>1.12*MCMC_GAL[k].BulgeMass[snap]-4.12)
        {
          if((MCMC_GAL[k].BlackHoleMass[snap]>-0.8929*MCMC_GAL[k].BulgeMass[snap]+15.8)
              && (MCMC_GAL[k].BlackHoleMass[snap]<-0.8929*MCMC_GAL[k].BulgeMass[snap]+17.88))
              binblackholeup[0]+=MCMC_GAL[k].Weight[snap];
          else if (MCMC_GAL[k].BlackHoleMass[snap]>-0.8929*MCMC_GAL[k].BulgeMass[snap]+17.88)
              binblackholeup[1]+=MCMC_GAL[k].Weight[snap];
        }
      else if((MCMC_GAL[k].BlackHoleMass[snap]>-0.8929*MCMC_GAL[k].BulgeMass[snap]+15.8)
              && (MCMC_GAL[k].BlackHoleMass[snap]<-0.8929*MCMC_GAL[k].BulgeMass[snap]+17.88))
              binblackholedown[0]+=MCMC_GAL[k].Weight[snap];
      else if(MCMC_GAL[k].BlackHoleMass[snap]>-0.8929*MCMC_GAL[k].BulgeMass[snap]+17.88)
              binblackholedown[1]+=MCMC_GAL[k].Weight[snap];
    }
  //printf("binup0=%f bindown0=%f binup1=%f bindown1=%f\n",
  //	  binblackholeup[0], binblackholedown[0], binblackholeup[1], binblackholedown[1]);
}


double chi_square_probability(int ObsNr, double *samdata, int snap)
{
  int df, i, knstrn = 0;
  double chsq = 0.0, temp = 0.0, prob = 0.0;
  double obs, obserror;

  df = Nbins[snap][ObsNr] - knstrn;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      obs = MCMC_Obs[ObsNr].Obs[snap][i];
      obserror = MCMC_Obs[ObsNr].Error[snap][i];

      if(obs < 0.0 || (obs == 0 && samdata[i] > 0))
    	printf("Bad expected number in chsone\n");
      if(obs == 0 && samdata[i] == 0)
	--df;
      else
        {
    	  if(snap==0)
    	    {
    		  if((obserror/obs)<MCMC_Minimum_Obs_Error)
    			obserror=MCMC_Minimum_Obs_Error*obs;
    	    }
    	  else
    	    {
    		  //if((obserror/obs)<3.*MCMC_Minimum_Obs_Error)
    			//obserror=3.*MCMC_Minimum_Obs_Error*obs;
    		  if((obserror/obs)< MCMC_Minimum_Obs_Error)
    		  	obserror= MCMC_Minimum_Obs_Error*obs;
    		  //if(ObsNr==1)
    		  //	obserror=5.*MCMC_Minimum_Obs_Error*obs;
    	    }

          temp=samdata[i]-obs;
          chsq += (temp*temp)/(obserror*obserror);
          //if(ObsNr==0)
          //printf("OBS[%d] snap=%d bin=%f sam=%f obs=%f error=%f chsq=%f\n",
          //		ObsNr, snap, MCMC_Obs[ObsNr].Bin_low[snap][i]+(MCMC_Obs[ObsNr].Bin_high[snap][i]-MCMC_Obs[ObsNr].Bin_low[snap][i])/2.,
          //		samdata[i], obs, obserror,(temp*temp)/(obserror*obserror));
        }
    }


  prob=exp(-chsq/2.);
  return prob;
}


double maximum_likelihood_probability(int ObsNr, double *samfract, int snap)
{
  double fracterr = 0.025, probaux = 99.0, prob = 99.0;
  int i;
  double obs, obserror;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      obs = MCMC_Obs[ObsNr].Obs[snap][i];
      obserror = MCMC_Obs[ObsNr].Error[snap][i];
      probaux=exp(-(pow2(samfract[i]-obs)/(2.*(fracterr*fracterr+obserror*obserror))));
      if(i==0) prob=probaux;
      else prob=prob*probaux;
    }

  return prob;
}

double binomial_probability(int ObsNr, double *samup, double *samdown, int snap)
{
  int i;
  double prob = 99.0, probaux = 99.0;
  double obsup, obsdown;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      obsup = MCMC_Obs[ObsNr].ObsUp[snap][i];
      obsdown = MCMC_Obs[ObsNr].ObsDown[snap][i];

      if((1.0-betai(obsup,obsdown+1.0,samup[i]/(samup[i]+samdown[i])))>betai(obsup,obsdown+1.0,samup[i]/(samup[i]+samdown[i])))
        probaux=2.0*betai(obsup,obsdown+1.0,samup[i]/(samup[i]+samdown[i]))+1e-20;
      else
        probaux=2.0*(1.0-betai(obsup,obsdown+1.0,samup[i]/(samup[i]+samdown[i])))+1e-20;

      if(i==0) prob=probaux;
      else prob=prob*probaux;
    }
  //printf("\nBinomial Probability=%f\n",prob);
  return prob;
}


double gammp(double a, double x)
{
  if(x < 0.0 || a <= 0.0)
    printf("bad args in gammp\n");
  if(x == 0.0)
    return 0.0;
  else if((int) a >= ASWITCH)
    return gammpapprox(a, x, 1);
  else if(x < a + 1.0)
    return gser(a, x);
  else
    return 1.0 - gcf(a, x);
}



double gammq(double a, double x)
{
  if(x < 0.0 || a <= 0.0)
    printf("bad args in gammq\n");
  if(x == 0.0)
    return 1.0;
  else if((int) a >= ASWITCH)
    return gammpapprox(a, x, 0);
  else if(x < a + 1.0)
    return 1.0 - gser(a, x);
  else
    return gcf(a, x);
}



double gser(double a, double x)
{
  double gln;
  double sum, del, ap;
  int n;

  gln = gammln(a);
  ap = a;
  del = sum = 1.0 / a;

  for(n = 1; n <= ASWITCH; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;
      if(fabs(del) < fabs(sum)*EPS) return sum*exp(-x+a*log(x)-gln);
    }
}




double gcf(double a, double x)
{
  double gln;
  int i;
  double an, b, c, d, del, h;

  gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;

  for(i = 1; i <= ASWITCH; i++)
    {
      an = -i * (i - a);
      b += 2.0;
      d = an * d + b;
      if(fabs(d) < FPMIN)
	    d = FPMIN;
      c = b + an / c;
      if(fabs(c) < FPMIN)
	    c = FPMIN;
      d = 1.0 / d;
      del = d * c;
      h *= del;
      if(fabs(del-1.0) <= EPS) break;
    }
  return exp(-x + a * log(x) - gln) * h;
}




double gammpapprox(double a, double x, int psig)
{
  double y[18] = { 0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
    0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
    0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
    0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
    0.87126389619061517, 0.95698180152629142
  };
  double w[18] = { 0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
    0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
    0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
    0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
    0.085346685739338721, 0.085983275670394821
  };
  int j, ngau = 18;
  double xu, t, sum, ans;
  double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
  double gln;

  gln = gammln(a);

  if(x > a1)
    if((a1 + 11.5 * sqrta1) > (x + 6.0 * sqrta1))
      xu = a1 + 11.5 * sqrta1;
    else
      xu = x + 6.0 * sqrta1;
  else if((a1 - 7.5 * sqrta1) < (x - 5.0 * sqrta1))
    if(0. > (a1 - 7.5 * sqrta1))
      xu = 0;
    else
      xu = (a1 - 7.5 * sqrta1);
  else if(0. > (x - 5.0 * sqrta1))
    xu = 0;
  else
    if(0.>( x - 5.0*sqrta1)) xu =0;
    else xu=( x - 5.0*sqrta1);


  sum = 0;
  for(j = 0; j < ngau; j++)
    {
      t = x + (xu - x) * y[j];
      sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
    }
  ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);

  if(psig)
    {
      if(ans > 0.0)
	{
	  return ans + 1.;
	}
      else
	return -ans;
    }
  else
    {
      if(ans >= 0.0)
	{
	  return ans;
	}
      else
	return ans + 1.;
    }
}




double gammln(double xx)
{
  int j;
  double x, tmp, y = 0, ser;

  double cof[14] = { 57.1562356658629235, -59.5979603554754912,
    14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
    .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
    -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
    .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5
  };
  if(xx <= 0)
    printf("bad arg in gammln\n");
  y = x = xx;
  tmp = x + 5.24218750000000000;
  tmp = (x + 0.5) * log(tmp) - tmp;
  ser = 0.999999999999997092;
  for(j = 0; j < 14; j++)
    ser += cof[j] / ++y;
  return tmp + log(2.5066282746310005 * ser / x);
}




double betai(double a, double b, double x)
{
  double bt;
  if(x < 0.0 || x > 1.0)
    printf("Bad x in routine betai\n");
  if(x == 0.0 || x == 1.0)
    bt = 0.0;
  else
    bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
  if(x < (a + 1.0) / (a + b + 2.0))
    return bt * betacf(a, b, x) / a;
  else
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}






double betacf(double a, double b, double x)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) printf("a or b too big, or MAXIT too small in betacf\n");
  return h;
}
#endif //MCMC

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef FPMIN
#undef ASWITCH
#undef MAXIT
