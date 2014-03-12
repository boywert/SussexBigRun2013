#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"


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
	double prob[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double final_prob, prob0;
	int i, snap;
	char buf[1000];
	FILE *fa;
 
	/* Bin samdata into binsamdata according to observational constraints.
	 * The first argument of the bin functions and of the likelihood
	 * function (Chi^2 ot MLM) indicates the observational data set to use.*/

	//FIRST CHI^2 tests then BINOMIAL!!!!

	//ALL Luminosity Function to be compared in units of M-5logh and phi(Mpc-3h-3mag-1)

	final_prob=1.;

	for(snap=0;snap<NOUT;snap++)
	  {
		//printf("totgals=%d\n",TotMCMCGals[snap]);
		if((samdata = malloc(sizeof(double) * TotMCMCGals[snap])) == NULL)
		  exit(999);

		if((binsamdata = malloc(sizeof(double) * Nbins[snap][0])) == NULL)
		  exit(999);
		for(i = 0; i < TotMCMCGals[snap]; i++)
	      samdata[i] = MCMC_GAL[i].StellarMass[snap];
		bin_function(0, binsamdata, samdata, snap);
		//stellar mass function
		prob[0]=chi_square_probability(0, binsamdata, snap);
		OutputObs[snap][0][CurrentMCMCStep][Nbins[snap][0]]=prob[0];
		free(binsamdata);

		if((binsamdata=malloc(sizeof(double)*Nbins[snap][1]))==NULL)
		  exit(999);
		for(i=0;i<TotMCMCGals[snap];i++)
		  samdata[i]=MCMC_GAL[i].MagK[snap];
		bin_function(1, binsamdata, samdata, snap);
		//K-band luminosity function
		prob[1]=chi_square_probability(1, binsamdata, snap);
		free(binsamdata);

		if((binsamdata=malloc(sizeof(double)*Nbins[snap][2]))==NULL)
		  exit(999);
		if((float)(int)(ZZ[ListOutputSnaps[snap]]+0.5)==0.) //Bj band at z=0
		  for(i=0;i<TotMCMCGals[snap];i++)
			  samdata[i]=MCMC_GAL[i].MagB[snap]-0.28*(MCMC_GAL[i].MagB[snap]-MCMC_GAL[i].MagV[snap]);
		else //BBand at all other z
		  for(i=0;i<TotMCMCGals[snap];i++)
			samdata[i]=MCMC_GAL[i].MagB[snap];
		bin_function(2, binsamdata, samdata, snap);
		//B-band luminosity function
		prob[2]=chi_square_probability(2, binsamdata, snap);
		free(binsamdata);

		if((binsamdata=malloc(sizeof(double)*Nbins[snap][8]))==NULL)
			exit(999);
		bin_colors(8, binsamdata, snap);
		prob[8]=maximum_likelihood_probability(8, binsamdata, snap);
		free(binsamdata);

	  bin_bhbm(binblackholeup, binblackholedown, snap);
	  prob[10]=binomial_probability(10, binblackholeup, binblackholedown, snap);

	  /*if((binsamdata=malloc(sizeof(double)*Nbins[snap][9]))==NULL) exit(999);
	  bin_color_hist(9, binsamdata);
	  prob[10]=chi_square_probability(10, binsamdata, snap);
	  free(binsamdata);*/

		/*if((binsamdata=malloc(sizeof(double)*Nbins[snap][8]))==NULL)
		  exit(999);
		  bin_colors(8, binsamdata, snap);
		  prob[8]=maximum_likelihood_probability(8, binsamdata, snap);
		  free(binsamdata);

		  if((binsamdata = malloc(sizeof(double) * Nbins[snap][9])) == NULL)
		  exit(999);
		  for(i = 0; i < TotMCMCGals[snap]; i++)
		  samdata[i] = MCMC_GAL[i].ColdGas[snap];
		  bin_function(9, binsamdata, samdata, snap);
		  prob[9]=chi_square_probability(9, binsamdata, snap);
	  	  free(binsamdata);



  	      if((binsamdata=malloc(sizeof(double)*Nbins[snap][2]))==NULL) exit(999);
  	      for(i=0;i<TotMCMCGals[snap];i++) samdata[i]=MCMC_GAL[i].Magu[snap];
  	      bin_function(2, binsamdata, samdata, snap);
  	      prob[2]=chi_square_probability(2, binsamdata, snap);
  	      free(binsamdata);

  	      if((binsamdata=malloc(sizeof(double)*Nbins[snap][3]))==NULL) exit(999);
  	      for(i=0;i<TotMCMCGals[snap];i++) samdata[i]=MCMC_GAL[i].Magg[snap];
  	      bin_function(3, binsamdata, samdata, snap);
  	      prob[3]=chi_square_probability(3, binsamdata, snap);
  	      OutputObs[snap][3][CurrentMCMCStep][Nbins[snap][3]]=prob[3];
  	      free(binsamdata);

  	      if((binsamdata=malloc(sizeof(double)*Nbins[snap][4]))==NULL) exit(999);
  	      for(i=0;i<TotMCMCGals[snap];i++) samdata[i]=MCMC_GAL[i].Magr[snap];
  	      bin_function(4, binsamdata, samdata, snap);
  	      prob[4]=chi_square_probability(4, binsamdata, snap);
  	      free(binsamdata);

  	      if((binsamdata=malloc(sizeof(double)*Nbins[snap][5]))==NULL) exit(999);
  	      for(i=0;i<TotMCMCGals[snap];i++) samdata[i]=MCMC_GAL[i].Magi[snap];
  	      bin_function(5, binsamdata, samdata, snap);
  	      prob[5]=chi_square_probability(5, binsamdata, snap);
  	      free(binsamdata);

  	      if((binsamdata=malloc(sizeof(double)*Nbins[snap][6]))==NULL) exit(999);
  	      for(i=0;i<TotMCMCGals[snap];i++) samdata[i]=MCMC_GAL[i].Magz[snap];
  	      bin_function(6, binsamdata, samdata, snap);
  	      prob[6]=chi_square_probability(6, binsamdata, snap);
  	      free(binsamdata);




  	      */

		free(samdata);

		printf("snap=%d\n prob[0]=%0.5e\n prob[1]=%0.5e\n prob[2]=%0.5e\n prob[8]=%0.5e\n prob[20]=%0.5e\n prob=%0.5e\n",
				snap, prob[0], prob[1], prob[2], prob[8], prob[10], prob[0]*prob[1]*prob[2]);

		if(snap==0)
			prob0=prob[0];
		/*if(snap==0)
		 final_prob*=prob[0]*prob[1]*prob[2]*prob[8]*prob[10];
		else if(snap<3)
			final_prob*=prob[0]*prob[1]*prob[2]*prob[8];
		else*/
			final_prob*=prob[0]*prob[1]*prob[2];

			//final_prob*=prob[0];
	  }
	printf("final prob=%0.5e\n",final_prob);

	if(final_prob>1.e-22 && prob0>1.e-4)
		printf("ola\n");

	return final_prob;
}


/**@brief Bin Luminosity or Stellar Mass Function*/
void bin_function(int ObsNr, double *binsamdata, double *samdata, int snap)
{
  int i, k;
  double Bin;
  FILE *fa;
  char buf[1000];

#ifndef PARALLEL
  sprintf(buf, "/galformod/scratch/bmh20/plots/data/mcmc_plus_obs%d_z%d.txt",ObsNr,(int)(ZZ[ListOutputSnaps[snap]]+0.5));
 	  if(!(fa = fopen(buf, "w")))
 	    {
 		  char sbuf[1000];
 		  sprintf(sbuf, "can't open file `%s'\n", buf);
 		  terminate(sbuf);
 	    }
#endif

  if((MCMC_Obs[1].Bin[snap][ObsNr] - MCMC_Obs[0].Bin[snap][ObsNr]) > 0.)
    Bin = MCMC_Obs[1].Bin[snap][ObsNr] - MCMC_Obs[0].Bin[snap][ObsNr];
  else
    Bin = MCMC_Obs[0].Bin[snap][ObsNr] - MCMC_Obs[1].Bin[snap][ObsNr];

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      binsamdata[i] = 0;

      for(k=0;k<TotMCMCGals[snap];k++)
    	  	if(samdata[k]>=(MCMC_Obs[i].Bin[snap][ObsNr]-Bin/2.)
    	  			&& samdata[k] <= (MCMC_Obs[i].Bin[snap][ObsNr]+Bin/2.)) binsamdata[i]+=MCMC_GAL[k].Weight[snap];
      binsamdata[i]/=(BoxSize*BoxSize*BoxSize*Bin);
      OutputObs[snap][ObsNr][CurrentMCMCStep][i]=binsamdata[i];
#ifndef PARALLEL
      fprintf(fa, "%g %g %g %g\n", MCMC_Obs[i].Bin[snap][ObsNr], MCMC_Obs[i].Obs[snap][ObsNr],
    		  	                   MCMC_Obs[i].Error[snap][ObsNr], binsamdata[i]);
#endif

      //printf("%g %g %g %g\n", MCMC_Obs[i].Bin[snap][ObsNr], MCMC_Obs[i].Obs[snap][ObsNr],
      // 	  	                   MCMC_Obs[i].Error[snap][ObsNr], binsamdata[i]);
    }
#ifndef PARALLEL
  fclose(fa);
#endif
}

/**@brief Bin fraction of red galaxies.
 *    The separation is done based on g-r color
 *    at z=0 and on SFR for higher z */
void bin_colors(int ObsNr, double *binredfraction, int snap)
{
  int i, k, IsRedGalaxy;
  double red, blue, Bin;
  double color;

  Bin = MCMC_Obs[1].Bin[snap][ObsNr] - MCMC_Obs[0].Bin[snap][ObsNr];

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      binredfraction[i] = 0.;
      red = 0.;
      blue = 0.;

      for(k=0;k<TotMCMCGals[snap];k++)
      	  {
    	    //color cut for z=0
    	    if((float)(int)(ZZ[ListOutputSnaps[snap]]+0.5)==0.0)
    	      {
    	    	color=MCMC_GAL[k].Magu[snap]-MCMC_GAL[k].Magr[snap];
    	    	IsRedGalaxy=(color > (2.06-0.244*tanh((MCMC_GAL[k].Magr[snap]+20.07)/1.09)));
    	      }
    	    //SFR cut for high redshift
    	    else
    	    	IsRedGalaxy=(MCMC_GAL[k].Sfr[snap]/pow(10.,MCMC_GAL[k].StellarMass[snap]) < 0.01);

    	    if(IsRedGalaxy == 1)
    	  		{
    	    		if(MCMC_GAL[k].StellarMass[snap]>=(MCMC_Obs[i].Bin[snap][ObsNr]-Bin/2.0)
    	    			&& MCMC_GAL[k].StellarMass[snap] <= (MCMC_Obs[i].Bin[snap][ObsNr]+Bin/2.0)) red+=MCMC_GAL[k].Weight[snap];
    	  		}
    	    else
    	    	{
    	    		if(MCMC_GAL[k].StellarMass[snap]>=(MCMC_Obs[i].Bin[snap][ObsNr]-Bin/2.0)
    	    			&& MCMC_GAL[k].StellarMass[snap] <= (MCMC_Obs[i].Bin[snap][ObsNr]+Bin/2.0)) blue+=MCMC_GAL[k].Weight[snap];
    	    	}
      	  }

      binredfraction[i] = red * 1.0 / (blue * 1.0 + red * 1.0);
      //printf("mass=%f obs =%f samcolors=%f\n",MCMC_Obs[i].Bin[snap][ObsNr], MCMC_Obs[i].Obs[snap][ObsNr], binredfraction[i]);
    }
}

/**@brief Bin SAM colours according to observations of Baldry2004*/
void bin_color_hist(int ObsNr, double *bincolorhist, int snap)
{
  int i, k;
  double Bin, color;



  Bin = MCMC_Obs[1].Bin[snap][ObsNr] - MCMC_Obs[0].Bin[snap][ObsNr];

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      bincolorhist[i] = 0;

      for(k=0;k<TotMCMCGals[snap];k++)
      {
    	  color=MCMC_GAL[k].Magg[snap]-MCMC_GAL[k].Magr[snap];
    	  if(color>=(MCMC_Obs[i].Bin[snap][ObsNr]-Bin/2.)
    			  && color <= (MCMC_Obs[i].Bin[snap][ObsNr]+Bin/2.)
    			  && MCMC_GAL[k].StellarMass[snap] > 9.0 && MCMC_GAL[k].StellarMass[snap] < 9.5)
    		  bincolorhist[i]=bincolorhist[i]+MCMC_GAL[k].Weight[snap];
      }
      //printf("%f %f %f\n",MCMC_Obs[i].Bin[snap][ObsNr], MCMC_Obs[i].Obs[snap][ObsNr],bincolorhist[i]);
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

  /*for(j=0;j<10;j++)
     printf("chimass=%f\n",obs[j]); */
  df = Nbins[snap][ObsNr] - knstrn;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      obs = MCMC_Obs[i].Obs[snap][ObsNr];
      obserror = MCMC_Obs[i].Error[snap][ObsNr];

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
    		  if((obserror/obs)<3.*MCMC_Minimum_Obs_Error)
    			obserror=3.*MCMC_Minimum_Obs_Error*obs;
    		  //if(ObsNr==1)
    		  //	obserror=5.*MCMC_Minimum_Obs_Error*obs;
    	    }

    	  //if(snap==3 && ObsNr==2.)
    		 // obserror=2.*MCMC_Minimum_Obs_Error*obs;

          temp=samdata[i]-obs;
          //printf("snap=%d sam=%f obs=%f error=%f\n",snap,samdata[i],obs, obserror);
          chsq += (temp*temp)/(obserror*obserror);
        }
    }

  prob = gammq(0.5 * (df), 0.5 * (chsq));
  return prob;
}


double maximum_likelihood_probability(int ObsNr, double *samfract, int snap)
{
  double fracterr = 0.1, probaux = 99.0, prob = 99.0;
  int i;
  double obs, obserror;

  for(i = 0; i < Nbins[snap][ObsNr]; i++)
    {
      obs = MCMC_Obs[i].Obs[snap][ObsNr];
      obserror = MCMC_Obs[i].Error[snap][ObsNr];
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
      obsup = MCMC_Obs[i].ObsUp[snap][ObsNr - NChiTests];
      obsdown = MCMC_Obs[i].ObsDown[snap][ObsNr - NChiTests];

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
