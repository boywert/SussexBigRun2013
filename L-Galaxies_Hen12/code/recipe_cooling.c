#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_cooling.c
 *  @brief recipe_cooling.c calculates the amount of mass that cools
 *  from the hot to the cold phase at each timestep (this fraction
 *  is then reduced by AGN heating) and updates the hot and cold gas
 *  fractions accordingly.
 *
 * This recipe calculates the amount of mass that cools from the hot
 * to the cold phase at each timestep. Two different infalling regimes
 * are assumed depending on the redshift and mass of the halos. This is
 * the recipe proposed by (White & Rees, 1978) and that has since been
 * used in most SA.
 *
 * Before only central galaxies could cool gas, so R_hot was just set
 * equal to Rvir. After Guo2010, satellites can have cooling, but since
 * they are loosing dark matter, Rvir is not a good approximation of R_Hot
 * so Gal.HotRadius was introduced.
 *
 * -> At early times and for low mass halos, the cooling radius can be
 * larger then the virial radius. In this case, despite the infalling
 * gas being shock heated to the virial temperature, it condenses within
 * a halo dynamical time and a quasi-static atmosphere cannot form.
 * Single line on the code, with no calculations needed, just that
 * \f$t_{\rm{dyn,h}}=R_{\rm{vir}}/V_{\rm{vir}}\f$
 * and therefore
 * \f$\dot{M}_{\rm{cool}}=0.5\frac{m_{\rm{hot}}}{t_{\rm{dyn}}}
 * =0.5m_{\rm{hot}}\frac{V_{\rm{vir}}}{R_{\rm{vir}}}\f$
 *
 *
 * -> For massive halos and at late times, the cooling radius lies within
 * the virial radius and the infalling gas does form a quasi-static hot
 * atmosphere that extends to the virial radius. This gas can cool at
 * later times and its accretion into central regions is modeled through
 * a cooling flow.
 *
 * \f$\dot{M}_{\rm{cool}}=
 *  m_{\rm{hot}}\frac{r_{\rm{cool}}}{R_{\rm{vir}}}\frac{1}{t_{\rm{cool}}}\f$
 *  eq5 (no 0.5 factor...)
 *
 *
 * In both cases the condensed gas settles into a central disk where
 * star formation occurs. The cooling rate is given by the cooling-AGN
 * heating.
 *
 * BH quiet accretion rate:
 *  \f$\dot{m}_{\rm BH,R}=k_{\rm AGN}\left(\frac{m_{\rm
    BH}}{10^8\,M_{\odot}}\right)\left(\frac{f_{\rm
    hot}}{0.1}\right)\left(\frac{V_{\rm
    vir}}{200\,\mathrm{km\,s}^{-1}}\right)^3\f$

   Luminosity from quiet accretion:
   \f$L_{\rm BH}=\eta\,\dot{m}_{\rm BH,R}\,c^2,\f$

   Corresponding reduction in cooling:
    \f$\dot{m}'_{\rm{cool}}=
    \dot{m}_{\rm{cool}}-\frac{L_{BH}}{\frac{1}{2}V^2_{\rm{vir}}}\f$
 *
 *
*/

/** @brief main cooling recipe, where the cooling rates are calculated */

double cooling_recipe(int p, double dt)
{
  double Vvir, Rvir, x, lambda, tcool, rcool, temp, tot_hotMass, tot_metals, HotRadius;
  double coolingGas, logZ, rho_rcool, rho0;

  mass_checks("cooling_recipe #1",p);

  tot_hotMass = Gal[p].HotGas;
#ifdef METALS
  //tot_metals = Gal[p].MetalsHotGas.type2;
  tot_metals = metals_total(Gal[p].MetalsHotGas);
#else
  tot_metals = Gal[p].MetalsHotGas;
#endif

  if(tot_hotMass > 1.0e-6) {
    /* TODO - Should Rvir be used at all in this recipe after Guo2010?
     * probably always HotRadius*/
    Vvir = Gal[p].Vvir;
    Rvir = Gal[p].Rvir;
    
      
      tcool = Rvir / Vvir; // tcool = t_dynamical = Rvir/Vvir


      /* temp -> Temperature of the Gas in Kelvin, obtained from
       * hidrostatic equilibrium KT=0.5*mu_p*(Vc)^2 assuming Vvir~Vc */
      temp = 35.9 * Vvir * Vvir;
      
    if (Gal[p].Type == 0)
      HotRadius = Gal[p].Rvir;
    else 
      HotRadius = Gal[p].HotRadius;
    
    if(tot_metals > 0)
      logZ = log10(tot_metals / tot_hotMass);
    else
      logZ = -10.0;

    //eq. 3 and 4 Guo2010
    lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
    x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
    x /= (UnitDensity_in_cgs * UnitTime_in_s);  // now in internal units
    rho_rcool = x / (0.28086 * tcool);
    /* an isothermal density profile for the hot gas is assumed here */
    rho0 = tot_hotMass / (4 * M_PI * HotRadius);
    rcool = sqrt(rho0 / rho_rcool);
    
    if (Gal[p].CoolingRadius < rcool)
      Gal[p].CoolingRadius = rcool;
      
    //if Hotradius is used, when galaxies become type 1's there will be a discontinuity in the cooling
    if(rcool > Rvir) // INFALL DOMINATED REGIME
      //coolingGas = tot_hotMass; - Delucia 2007
      /*comes in to keep the continuity (Delucia2004) */
      //put h
      coolingGas = tot_hotMass / (HotRadius / Vvir) * dt;
    else // HOT PHASE REGIME -> TODO - WHERE IS THE 0.5 factor of eq 5
      /*coolingGas = (tot_hotMass / Rvir) * (rcool / tcool) * dt * 0.5; */
      coolingGas = (tot_hotMass / HotRadius) * (rcool / tcool) * dt ;
        
    //Photoionizing background
    if (log10(temp) < 4.0)
      coolingGas = 0.;
    
    if(coolingGas > tot_hotMass)
      coolingGas = tot_hotMass;
    else if(coolingGas < 0.0)
      coolingGas = 0.0;      
    
    /*  Suppress cooling due to AGN feedback, the gas is not actually heated,
     *  just the amount of cooling reduced. */
  mass_checks("cooling_recipe #1.5",p);
    if(AGNrecipeOn > 0)
      coolingGas -= do_AGN_heating(coolingGas, p, dt, x);
  mass_checks("cooling_recipe #2",p);

    if(coolingGas < 0.0)
      coolingGas = 0.0;  
  }
  else {
    coolingGas = 0.0;
  }
  
  /* determine the xray luminosity of any cooling gas in this snapshot
   * (White & Frenk 1991 eq21) */
  if(coolingGas > 0.0)
    Gal[p].XrayLum =
      log10(2.5 * (coolingGas / dt) * 6.31 * Gal[p].Vvir * Gal[p].Vvir) + 35.0;
  else
    Gal[p].XrayLum = 0.0;
  Gal[p].CoolingGas += coolingGas;
  
  return coolingGas;

}   



/** @brief calculates the energy released by black holes due to passive accretion,
  * that will be used to reduced the cooling.*/

double do_AGN_heating(double coolingGas, int centralgal, double dt, double x)
{
  double AGNrate, AGNheating, AGNaccreted, AGNcoeff, fraction, EDDrate, FreeFallRadius;

  /** @brief do_AGN_heating calculates the amount of energy released by
    * black holes due to passive accretion, Which is then used to reduced
    * the cooling.
    *
    * There is one parameter, AgnEfficiency, which is the efficiency of AGN
    * passive accretion and consequently of cooling flow reheating, Eqs. 10,
    * 11 & 12 in Croton 2006. The standard recipe is  AGNrecipeOn =1 (empirical)
    * with options 2 and 3 representing Bondi-Hoyle and cold cloud accretion.
    * The three should be identical and the use of empirical avoids people
    * shouting about duty-cycles beying inconsistent with Bondi-Hoyle & etc.
    */

  if(Gal[centralgal].HotGas > 0.0)
    {

      if(AGNrecipeOn == 2)
        {
          /* Bondi-Hoyle accretion recipe -- efficiency = 0.15
           * Eq. 29 in Croton 2006 */
          AGNrate = (2.5 * M_PI * G) * (0.75 * 0.6 * x) * Gal[centralgal].BlackHoleMass * 0.15;
        }
      else if(AGNrecipeOn == 3)
        {
          /* Cold cloud accretion recipe -- trigger: Rff = 50 Rdisk,
           * and accretion rate = 0.01% cooling rate
           * Eq. 25 in Croton 2006 */
          FreeFallRadius =
            Gal[centralgal].HotGas / (6.0 * 0.6 * x * Gal[centralgal].Rvir * Gal[centralgal].Vvir) / Gal[centralgal].HotRadius * Gal[centralgal].Rvir;
          if(Gal[centralgal].BlackHoleMass > 0.0 && FreeFallRadius < Gal[centralgal].GasDiskRadius * 50.0)
            AGNrate = 0.0001 * coolingGas / dt;
          else
            AGNrate = 0.0;
        }
      else
        {
            //empirical (standard) accretion recipe - Eq. 10 in Croton 2006 
            AGNrate =
            AgnEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
            * (Gal[centralgal].BlackHoleMass / 0.01) * pow3(Gal[centralgal].Vvir / 200.0)
            * ((Gal[centralgal].HotGas / Gal[centralgal].HotRadius * Gal[centralgal].Rvir / Gal[centralgal].Mvir) / 0.1);

         //optimized
	  /*   AGNrate = AgnEfficiency * (SOLAR_MASS * UnitTime_in_s / UnitMass_in_g * SEC_PER_YEAR)
	   * Gal[centralgal].BlackHoleMass * 100. * pow3(Gal[centralgal].Vvir * 0.005)
	   * Gal[centralgal].HotGas * Gal[centralgal].Rvir * 10. / (Gal[centralgal].Mvir * Gal[centralgal].HotRadius);*/

        }

      /* Eddington rate */
      /* Note that this assumes an efficiency of 50% 
       * - it ignores the e/(1-e) factor in L = e/(1-e) Mdot c^2 */
      EDDrate = 1.3e48 * Gal[centralgal].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;

      /* accretion onto BH is always limited by the Eddington rate */
      if(AGNrate > EDDrate)
        AGNrate = EDDrate;

      /*  accreted mass onto black hole */
      AGNaccreted = AGNrate * dt;

      /* cannot accrete more mass than is available! */
      if(AGNaccreted > Gal[centralgal].HotGas)
        AGNaccreted = Gal[centralgal].HotGas;

      /*  coefficient to heat the cooling gas back to the virial temperature of the halo */
      /*  1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s
       *  Eqs. 11 & 12 in Croton 2006 */
      AGNcoeff = (1.34e5 / Gal[centralgal].Vvir) * (1.34e5 / Gal[centralgal].Vvir);

      /*  cooling mass that can be suppressed from AGN heating */
      AGNheating = AGNcoeff * AGNaccreted;

      /* limit heating to cooling rate */
      if(AGNheating > coolingGas)
        {
          AGNaccreted = coolingGas / AGNcoeff;
          AGNheating = coolingGas;
        }

      /*  accreted mass onto black hole */
      fraction=AGNaccreted/Gal[centralgal].HotGas;
      Gal[centralgal].BlackHoleMass += AGNaccreted;
      Gal[centralgal].RadioAccretionRate += AGNaccreted / (dt*STEPS);
      Gal[centralgal].HotGas -= AGNaccreted;
      Gal[centralgal].MetalsHotGas = 
	metals_add(Gal[centralgal].MetalsHotGas,Gal[centralgal].MetalsHotGas,
		   -fraction);
#ifdef YIELDS
      Gal[centralgal].HotGas_elements = elements_add(Gal[centralgal].HotGas_elements,Gal[centralgal].HotGas_elements,-fraction);
#endif
#ifdef METALS_SELF
      Gal[centralgal].MetalsHotGasSelf = 
	metals_add(Gal[centralgal].MetalsHotGasSelf,Gal[centralgal].MetalsHotGasSelf,
		   -fraction);
#endif	
    }

  else
    AGNheating = 0.0;

  return AGNheating;

}


/** @brief updates the fractions of hot and cold gas due to cooling. */

void cool_gas_onto_galaxy(int p, double mcool)
{
  double fraction,Mdisk;
  int i;

  /** @brief cool_gas_onto_galaxy updates the fractions of hot and cold gas
    * due to cooling. This is done for the mass, metals and, after Guo2010,
    * spin components */

  Mdisk=Gal[p].ColdGas;
  if (mcool>Gal[p].HotGas)
    mcool = Gal[p].HotGas;

  /*  add the fraction 1/STEPS of the total cooling gas to the cold disk */
  if(mcool > 0.0) { 
     
    // We already know that 0<mcool<=Gal[p].HotGas
    fraction=mcool/Gal[p].HotGas;
    transfer_gas(p,"Cold",p,"Hot",fraction);

    if (DiskRadiusMethod == 2) {
      if (Gal[p].ColdGas != 0.0)
	for (i=0;i<3;i++)
	  Gal[p].GasSpin[i]=(Gal[p].GasSpin[i]*Mdisk+Gal[p].HaloSpin[i]*mcool)/(Gal[p].ColdGas);
	  get_gas_disk_radius(p);
    }
  }
  
}
