#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_disrupt.c
 *  @brief recipe_disrupt.c checks if a type 2 satellite galaxy should
 *         or not be disrupted due to tidal forces
 *
 *  This routine takes into account the tidal effects that satellite
 *  galaxies experience while orbiting a central companion. Since the
 *  baryonic component is more compact and denser that the dark matter
 *  it assumes that only type 2 satellites are affected (those that have
 *  already lost their dark matter halo). It assumes that the disruption
 *  is complete and instantaneous with all the satellite material being
 *  transferred into the central galaxy.
 *
 *  The satellite is assumed to orbit a singular isothermal potential:
 *
 *  \f$\phi(R)=V^2_{\rm{vir}}\rm{ln} R\f$ (Eq. 28 Guo2010)
 *
 *  Assuming conservation of energy and angular momentum along the
 *  orbit, its pericentric distance (closest point to the centre along
 *  the orbit) can be estimated from:
 *
 *  \f$\left(\frac{R}{R_{\rm{peri}}}\right)^2=
 *  \frac{lnR/R_{\rm{peri}}+\frac{1}{2}(V/V_{\rm{vir}})^2}
 *  {\frac{1}{2}(V_{\rm{t}}/V_{\rm{vir}})^2}\f$
 *  (Eq. 29 Guo2010).
 *
 *  The main halo density at this point is compared with the baryonic
 *  mass (cold gas + stellar) density of the satellite within its half
 *  mass radius. If
 *
 *  \f$ \frac{M_{\rm{DM,halo}}(R_{\rm{peri}})}{R^3_{\rm{peri}}}\equiv
 *  \rho_{\rm{DM,halo}}>
 *  \rho_{\rm{sat}}\equiv\frac{M_{\rm{sat}}}{R^3_{\rm{sat,half}}}\f$
 *  (Eq. 30 Guo2010)
 *
 *  the galaxy is disrupted.
 *
 *  TODO: shouldn't galaxy p be voided after disruption?
 *
 *  */


void disrupt(int p, int centralgal)
{
  double rho_sat, rho_cen;
  double cen_mass, r_sat, radius;
  int outputbin, j, q;

  /* If the main halo density at the pericentre (closest point in the orbit
   * to the central galaxy)is larger than the satellite's density at the
   * half mass radius, the satellite is completely disrupted. Note that since
   * the satellite is a type 2 the only mass components remaining and
   * contributing to the density are the cold gas and stellar mass. */

  //TODO If we are passing in centralgal then we should not set it here 
  centralgal=Gal[p].CentralGal;
 
  mass_checks("Top of disrupt",centralgal);
  mass_checks("Top of disrupt",p);

  /* Radius calculated at the peri-point */
  radius=peri_radius(p, centralgal);
  if (radius < 0) {
   terminate("must be wrong \n");
  }

  /* Calculate the density of the main central halo at radius (the peri-centre).
   * The tidal forces are caused by the dark matter of the main halo, hence Mvir
   * is used. */
  cen_mass=Gal[centralgal].Mvir*radius/Gal[centralgal].Rvir;
  rho_cen=cen_mass/pow3(radius);

  /* Calculate the density of the satellite's baryonic material */
  if (Gal[p].DiskMass+Gal[p].BulgeMass>0)
  {
    /* Calculate the rho according to the real geometry */
    r_sat = sat_radius(p);

    /* Or use radius at the mean radius of the stellar mass */
    /* r_sat=(Gal[p].BulgeMass*Gal[p].BulgeSize+(Gal[p].StellarMass-Gal[p].BulgeMass)
     *        *Gal[p].StellarDiskRadius/3*1.68)
     *       /(Gal[p].StellarMass); */

    /* to calculate the density */
    //rho_sat=Gal[p].StellarMass/pow2(r_sat);
    rho_sat=(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas)/pow3(r_sat);
  }
  else
    rho_sat=0.0;

  /* If density of the main halo is larger than that of the satellite baryonic
   * component, complete and instantaneous disruption is assumed. Galaxy becomes
   * a type 3 and all its material is transferred to the central galaxy. */
  if (rho_cen > rho_sat) {
    Gal[p].Type = 3;
#ifdef GALAXYTREE      
    q = Gal[Gal[p].CentralGal].FirstProgGal;
    if (q >= 0)
    {
      // add progenitors of Gal[p] to the list of progentitors of Gal[p].CentralGal
      while (GalTree[q].NextProgGal >= 0)
      	q = GalTree[q].NextProgGal;
	
      GalTree[q].NextProgGal = Gal[p].FirstProgGal;
      
      if(GalTree[q].NextProgGal >= NGalTree)
	    {
	      printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		     q, p, GalTree[q].NextProgGal, NGalTree);
	      terminate("problem");
	    }
    }

    if(q < 0)
    	terminate("this shouldn't happen");
	
     // TODO if !(q>=0) shouldn't we set
      // Gal[Gal[p].CentralGal].FirstProgGal to Gal[p].FirstProgGal  ??

    q = GalTree[q].NextProgGal;

    if(q < 0)
    	terminate("inconsistency");

    if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
    	terminate("inconsistency");

    HaloGal[GalTree[q].HaloGalIndex].DisruptOn = 1;
#endif
    /* Put gas component to the central galaxy hot gas and stellar material into the ICM.
     * Note that the staellite should have no extended components. */
    // TODO Shouldn't the stars end up in the bulge? - this is a close merger
    transfer_gas(centralgal,"Hot",p,"Cold",1.);
    transfer_gas(centralgal,"Hot",p,"Hot",1.);
#ifdef TRACK_BURST
    /* Transfer burst component first */
    transfer_stars(centralgal,"Burst",p,"Burst",
		   (Gal[p].DiskMass+Gal[p].BulgeMass)/(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ICM));
#endif
    transfer_stars(centralgal,"ICM",p,"Disk",1.);
    transfer_stars(centralgal,"ICM",p,"Bulge",1.);
    /* Add satellite's luminosity into the luminosity of the ICL
     * component of the central galaxy. */
#ifdef ICL
    for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++)
      {
#ifdef OUTPUT_REST_MAGS 
      	Gal[centralgal].ICLLum[j][outputbin] += Gal[p].Lum[j][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
      	Gal[centralgal].ObsICL[j][outputbin] += Gal[p].ObsLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
      	Gal[centralgal].dObsICL[j][outputbin] += Gal[p].dObsLum[j][outputbin];
#endif
#endif
      }  
    }
#endif      
  } //if (rho_cen > rho_sat)
  mass_checks("Bottom of disrupt",centralgal);
  mass_checks("Bottom of disrupt",p);
  
}

/** @brief Calculates the distance of the satellite to the pericentre of the
  *        main dark matter halo. */

double peri_radius(int p, int centralgal)
{
  int i;
  double a, b, v[3], r[3], x, x0;
  for(i = 0; i < 3; i++)
    {
	  r[i] = wrap(Gal[p].Pos[i]-Gal[centralgal].Pos[i],BoxSize);
	  r[i] /= (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	  v[i] = Gal[p].Vel[i] - Gal[centralgal].Vel[i];
    }

  b = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / pow2(Gal[centralgal].Vvir);
  a =
    1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] -
	      pow2(r[0] * v[0] + r[1] * v[1] + r[2] * v[2])
		   / (r[0] * r[0] + r[1] * r[1] + r[2] * r[2])) / pow2(Gal[centralgal].Vvir);

  x = sqrt(b / a);
  x0 = 1000;
  while(abs(x0 - x) >= 1.e-8)
    {
      x0 = x;
      x = sqrt((log(x0) + b) / a);
    }
  if(x == 0)
    {
      terminate("wrong in peri_radius \n");
    }

  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/x;
}


/** @brief Calculates the half mass radius of satellite galaxies */
double sat_radius(int p)
{
  double r, rd, rb, Mdisk, Sigma0, rmin, rmax, rbin, M, tmprmax;
  double Mgas, Mbulge, rgd, Sigma0_g, rmi, rma, dr, totmass;
  int N = 50, i;
  //disk profile, rd is known
  rd=Gal[p].StellarDiskRadius/3.;
  rb=Gal[p].BulgeSize;
  Mdisk=Gal[p].DiskMass;
  rgd = Gal[p].GasDiskRadius/3.;
  Mgas = Gal[p].ColdGas;
  Mbulge = Gal[p].BulgeMass;
  totmass = Gal[p].DiskMass + Gal[p].BulgeMass + Gal[p].ColdGas;


  //defining a maximum search radius for half mass?
  tmprmax=(rd*1.68 > rb)? rd*1.68:rb;
  rmax = (tmprmax > rgd * 1.68)? tmprmax : rgd * 1.68;
  
  i = 0;
  rmin = 1.e-6;
  if (rmax == 0 || rmax < rmin)
    return(rmax);
  rbin=(rmax-rmin)/N;
  
  if (rd == 0 && rgd == 0)
    return(rb);
  if (rd ==0 && rb == 0)
    return(rgd*1.68);
  if (rb==0 && rgd==0)
    return(rd*1.68);
  
  /* Disk surface densities, used to calculate the total mass inside
   * a certain radius. */
  Sigma0 = Mdisk / (2 * M_PI * rd *rd);
  Sigma0_g = Mgas / (2 * M_PI *rgd * rgd);
  if (rgd == 0) {
    Sigma0_g =0;
    rgd = 1.;
  }
  if (rd == 0) {
    Sigma0 = 0;
    rd = 1.;
  }

  /* increases the search radius until it encompasses half the total mass taking
   * into account the stellar disk, stellar bulge and cold gas disk. */
  M = 0.0;
  do {
    rmi = (rmin) + i* rbin;
    rma = (rmin) + rbin * (i+1);
    dr = rma - rmi;
    r = (rma + rmi) /2; 
    M = M + diskmass_r(r,rd,Sigma0,dr)+bulgemass_r(r,rb,Mbulge,dr)+diskmass_r(r,rgd,Sigma0_g,dr);
    i++;
  }
  while(M < 0.5*totmass);
  return (r);
}


/** @brief Returns the mass of a disk at a certain radius.
 *         Disk profile -> exponential */
double diskmass_r(double r, double rd, double Sigma0,double dr)
{
  return 2 * M_PI * Sigma0 * dr * r * exp(-r / rd);
}


/** @brief Returns the mass of a bulge at a certain radius.
 *         Bulge profile -> de Vaucouleurs type r^{1/4} law */
double bulgemass_r(double r, double rb, double Mbulge,double dr)
{
  return  Mbulge / (4 * M_PI * pow3(rb)) * 1 / pow2(r / rb) * 1 / pow2(1 + r/rb) * 4 * M_PI * r * r *dr;
}
