#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


/**@file recipe_stripping.c
 * @brief deal_with_satellites.c 
 *
 *       This is where the gas and ICM components of newly accreted satellites
 *       are treated.
 *
 *       There are basically 2 options for the way satellite components are
 *       added into centrals:
 *
 *         SatelliteRecipe ==0 
 *         If inside Rvir, hot and ejected gas from satellites of both types 1 and 2
 *         is instantaneously striped and added to type 0.
 *
 *         SatelliteRecipe ==1
 *         Type 1's keep an ejected component.
 *         Type 1's are stripped of hot and ejected gas gradually and later in the code.
 *         A fraction of the hot and ejected gas in the type 2's is
 *         added to the type 1 and the rest to the type 0.
 *         If satellites are outside Rvir, type 1 keeps all its components and receives
 *         everything from type 2's.
 *
 *         In these routines, centralgal is the type 0 at the centre of the halo;
 *         Gal[i].CentralGal is the galaxy around which each satellite orbits;
 *         For simplicity of reference in the comments in the code below, the latter
 *         will be called the type 1, even though it may be the same galaxy as the type 0.
 *
 *         */

void deal_with_satellites(int centralgal, int ngal)
{
  int i;
  double dis, fraction;

  for(i = 0; i < ngal; i++)    /* Loop over all galaxies in the FoF-halo */
    {
	  mass_checks("Top of deal_with_satellites a",i);
	  mass_checks("Top of deal_with_satellites b",centralgal);
	  mass_checks("Top of deal_with_satellites c",Gal[i].CentralGal);

	  /* dis is the separation of the type 1 from the type 0 */
	  dis=separation_gal(centralgal,Gal[i].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	
	  /* Instantaneous stripping of gas from satellites and no ejection of type 2 into type 1,
	   * still there is the condition on Rvir that determines that if a galaxy is a newly
	   * accreted type 2 outside Rvir of type 0, its gas will go into the type 1. If it's
	   * a type 1 outside Rvir of type 0, it retains all its gas.*/
	  if(SatelliteRecipe == 0)
        {
		  /* If galaxy is a satellite inside Rvir it will lose its hot and
		   * ejected gas into the hot gas component of the centralgal.
		   * Only galaxies within Rvir contribute to the central halo.*/
		  if ( dis < Gal[centralgal].Rvir && i != centralgal)
    	    {
			  transfer_gas(centralgal,"Hot",i,"Hot",1.);
			  /* Note: the original code transferrred from ejected to hot here, but had all other
			   * such transfers between equivalent phases.  I have altered for consistency. */
			  transfer_gas(centralgal,"Ejected",i,"Ejected",1.);
			  transfer_stars(centralgal,"ICM",i,"ICM",1.);
#ifdef ICL
			  transfer_ICL(centralgal, i, 1.);
#endif
			  Gal[i].HotRadius =0.;
    	    }
		  /* If its a type 1 outside Rvir it retains all its gas components, so do nothing
		   * else if (Gal[i].Type ==1) {}
		   * If galaxy is a type 2 outside Rvir of type 0, then all its gas components
		   * will be added to the type 1. */
		  else
			if (Gal[i].Type == 2)
			  {
				transfer_gas(Gal[i].CentralGal,"Hot",i,"Hot",1.);
				transfer_gas(Gal[i].CentralGal,"Ejected",i,"Ejected",1.);
				transfer_stars(Gal[i].CentralGal,"ICM",i,"ICM",1.);
#ifdef ICL
				transfer_ICL(Gal[i].CentralGal,i,1.);
#endif
				Gal[i].HotRadius =0.;
			  }
        }

	  /* SatelliteRecipe == 1 => Guo2010 non instantaneous treatment of
	   * gas stripping from satellites (this means ejected and hot gas of type 2's can
	   * be split between type 0 and type 1)*/
	  else if (SatelliteRecipe == 1)
        {
		  /* Stripping of Type 2 galaxies */
		  if (Gal[i].Type ==2)
    	    {
			  /* Only galaxies within rvir contribute to total baryon mass in central halo */
			  if (dis < Gal[centralgal].Rvir)
				/* The fraction of the lost material that goes to the immediate host
				 * (for abbreviation below call this the type 1 and the galaxy at the centre
    			 * of the main halo the type 0; however these could be the same galaxy.) */
    			fraction=Gal[Gal[i].CentralGal].HotRadius / Gal[Gal[i].CentralGal].Rvir;
			  else
    			fraction=1.;

			  Gal[i].HotRadius = 0.0;
			  /* Split type 2 hot and ejected gas between type 0 and type 1. */
			  // A fraction of the gas goes to the type 1
			  transfer_gas(Gal[i].CentralGal,"Hot",i,"Hot",fraction);
			  transfer_gas(Gal[i].CentralGal,"Ejected",i,"Ejected",fraction);
			  mass_checks("deal_with_satellites i #0",i);
			  mass_checks("deal_with_satellites Gal[i].CentraGal #0",Gal[i].CentralGal);
			  transfer_stars(Gal[i].CentralGal,"ICM",i,"ICM",fraction);
			  mass_checks("deal_with_satellites i #1",i);
			  mass_checks("deal_with_satellites Gal[i].CentraGal #1",Gal[i].CentralGal);
#ifdef ICL
			  transfer_ICL(Gal[i].CentralGal,i,fraction);
#endif  
			  // All the rest goes to the type 0
			  if (fraction < 1.)
    		    {
				  transfer_gas(centralgal,"Hot",i,"Hot",1.);
				  transfer_gas(centralgal,"Ejected",i,"Ejected",1.);
				  transfer_stars(centralgal,"ICM",i,"ICM",1.);
				  mass_checks("deal_with_satellites #2",i);
				  mass_checks("deal_with_satellites #2",centralgal);
#ifdef ICL
				  transfer_ICL(centralgal,i,1.);
#endif  
    		    }
    	    }
		  else if (Gal[i].Type == 1 && dis < Gal[centralgal].Rvir && Gal[i].HotGas > 0.0)
    	    {
			  /* If galaxy is a type 1 and its already inside Rvir of type 0, then the
			   * gradual stripping of the gas is computed */
			  /* Fraction of hotgas that the type 1 will lose.
			   * Note that hot_retain_sat also re-evaluates HotRadius*/
			  fraction=1.-hot_retain_sat(i,centralgal)/Gal[i].HotGas;
			  if (fraction < 0.)
    		    {
				  printf("***Error in hot_retain_sat - returns value larger than HotGas***\n");
				  exit(1);
    		    }

			  transfer_gas(centralgal,"Hot",i,"Hot",fraction);
			  transfer_gas(centralgal,"Ejected",i,"Ejected",fraction);
			  mass_checks("deal_with_satellites #3",i);
			  mass_checks("deal_with_satellites #3",centralgal);
			  transfer_stars(centralgal,"ICM",i,"ICM",fraction);
			  mass_checks("deal_with_satellites #4",i);
			  mass_checks("deal_with_satellites #4",centralgal);
#ifdef ICL
			  transfer_ICL(centralgal,i,fraction);
#endif
    	    }

		  mass_checks("deal_with_satellites #5",i);
		  mass_checks("deal_with_satellites #5",centralgal);
        }//end of SatelliteRecipe == 1


	  mass_checks("Bottom of deal_with_satellites i",i);
	  mass_checks("Bottom of deal_with_satellites centralgal",centralgal);
    }
  /* End of choice of Satellite recipes */

   return;
  
}


/** Gradual stripping of hot and ejected gas from type 1 satellites. 
 *  This is caused both by tidal and ram-pressure stripping.
 *  This function returns the actual mass of hot gas that the
 *  type 1 retains.
 *
 *  TIDAL STRIPPING
 *  Hot gas is tidally stripped at the same rate at which dark matter is
 *  stripped:
 *
 * \f$ \frac{M_{\rm{hot}}(R_{\rm{tidal}})}{M_{\rm{hot,infall}}}=
 *  \frac{M_{\rm{DM}}}{M_{\rm{DM,infall}}}\f$
 *
 *  Since the hot gas distribution is assumed to be \f$ \rho \propto r^{-2}\f$
 *  this means \f$ M_{\rm{hot}}(r) \propto r.\f$ Therefore, the tidal
 *  radius beyond gas is stripped is given by:
 *
 *  \f$ R_{\rm{tidal}}=
 *  \left(\frac{M_{\rm{DM}}}{M_{\rm{DM,infall}}}\right)R_{\rm{DM,infall}}\f$
 *
 *  RAM PRESSURE STRIPING
 *  Let \f$R_{r.p.}\f$ represent the distance from the centre of the satellite
 *  at which ram pressure striping equals its self-gravity. Then:
 *
 *  \f$ \rho_{\rm{sat}}(R_{\rm{r.p.}})V^2_{\rm{sat}}=
 *      \rho_{\rm{par}}(R_{\rm{orbit}})V^2_{\rm{orbit}}\f$
 *  Where the four terms represent respectively the density of the satellite
 *  at \f$R_{\rm{r.p.}}\f$, the virial velocity of the satellite at infall,
 *  the density of the parent halo at the radius of the satellite and the
 *  orbit velocity of the satellite (given by \f$V_{\rm{c}} of the parent halo\f$)
 *
 *  The stripping radius is given by
 *
 *  \f$R_{\rm{strip}}=min(R_{\rm{tidal}},R_{\rm{r.p.}})\f$
 *
 * */
double hot_retain_sat(int i, int centralgal)
{
  double hotremain;
  double R_2,R_1,Rt;
  double RetainFrac;
  
  if (Gal[centralgal].Type != 0) exit(0);

  /*Calculate tidal stripping radius*/
   Rt=Gal[i].Len*PartMass/Gal[i].Mvir*Gal[i].Rvir;

  /*Ram pressure stripping radius calculation*/
  /*First calculate the orbital radius of the satellite R_orbit*/
  /*R_1=sqrt(pow(Gal[centralgal].Pos[0]-Gal[i].Pos[0],2.)+
	   pow(Gal[centralgal].Pos[1]-Gal[i].Pos[1],2.)+
	   pow(Gal[centralgal].Pos[2]-Gal[i].Pos[2],2.))
    /(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);*/

   R_1=separation_gal(centralgal,i)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

  /*If the central galaxy has no hot gas, it exerts no ram pressure stripping on the
   * satellite. */
  if (Gal[centralgal].HotGas<1.e-6)
    R_2=Gal[i].HotRadius;
  else {
    /*R_r.p./R_orbit*/
    RetainFrac=sqrt(Gal[i].HotGas*Gal[centralgal].Mvir/(Gal[i].Mvir*Gal[centralgal].HotGas))
      *(Gal[i].Mvir/Gal[i].Rvir)/(Gal[centralgal].Mvir/Gal[centralgal].Rvir);
    /*R_r.p. */
    R_2=RetainFrac*R_1;
  }

  /*Get the smaller of tidal and ram pressure stripping radii.*/
  if (R_2 > Rt) R_2=Rt;
  
  /*if the stripping radius is larger then hot radius there is
   * no stripping*/
  if (R_2>Gal[i].HotRadius || Gal[i].HotGas < 1.e-8)
    hotremain=Gal[i].HotGas;	 
  /* If stripping radius is smaller than the hot radius */
  else { 
    /* Assuming M_hot(r) proportional to r, the remaining hot gas
     * is given by: */
    hotremain=Gal[i].HotGas*R_2/Gal[i].HotRadius;
    /* hot radius is updated to the stripping radius*/
    Gal[i].HotRadius=R_2;
    /* Check that HotRadius has sensible values */
    if (Gal[i].HotRadius < 1.e-8)
      Gal[i].HotRadius = Gal[i].Len*PartMass/Gal[i].Mvir*Gal[i].Rvir;
    if (Gal[i].HotRadius > Gal[i].Rvir)	  
      Gal[i].HotRadius = Gal[i].Rvir;
  }

  return hotremain;
}
