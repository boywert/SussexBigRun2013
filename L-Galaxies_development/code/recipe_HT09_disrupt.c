#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <gsl/gsl_rng.h>
//#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

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


double get_merging_radius(int halonr, int mother_halonr)
{
  int central_halonr;
  double SatelliteRadius, MotherHaloRvir;

  /*  recipe updated for more accurate merging time (see BT eq 7.26),
     now satellite radius at previous timestep is included */

  central_halonr = Halo[Halo[halonr].Descendant].FirstProgenitor;
  if(central_halonr == halonr)
    {
      printf("can't be...!\n");
      exit(1);
    }

  /* in physical length  (/(1 + ZZ[Halo[halonr].SnapNum]))*/
  SatelliteRadius = separation_halo(central_halonr,halonr)/(1 + ZZ[Halo[halonr].SnapNum]);

  MotherHaloRvir = get_virial_radius(mother_halonr);
  if(SatelliteRadius > MotherHaloRvir)
    SatelliteRadius = MotherHaloRvir;

  return SatelliteRadius;
}





double get_deltar(int galID, double deltaT)
{
  double coulomb, deltar, SatelliteMass, Bx, CentralMass, x, CentralVelDisp;

  //G is the Gravitational constant in Mpc Solar Masses^-1 (km/s)^2 in units of 10^10

  SatelliteMass=Gal[galID].Mergmass+Gal[galID].DiskMass+Gal[galID].BulgeMass+Gal[galID].ColdGas+Gal[galID].HotGas;
  CentralVelDisp=Gal[galID].CentralVelDisp;
  CentralMass=2.0*CentralVelDisp*CentralVelDisp*Gal[galID].MergRadiusconst/(G);
  x=1.0;
  //Bx=erff(x)-(2.0*x/sqrt(M_PI))*exp(-(x*x));
  Bx=1.;
  coulomb=1.0+CentralMass/SatelliteMass;

  //deltar gives the decay in the orbit of the satellite. MergerTimeMultiplier is defined to be
  //an increase in merger time if positive, so it should decrease deltar
  if(SatelliteMass > 0.0)
    {
      deltar=1./MergerTimeMultiplier*(G*SatelliteMass*log(coulomb)*Bx)*deltaT/
              (sqrt(2.0)*Gal[galID].CentralVelDisp*Gal[galID].MergRadius);
    }
  else
    deltar =0.0;

 return deltar;

}

void disruption_code (int galID, double time)
{
  float DisruptionRadius, CentralMass, CentralVelDisp;
  float Rd, a, Re, LostFract, BulgeVelDisp;
  float InitialBulgeMass, InitialDiskMass, InitialStellarMass;
  int j;

  //Values before disruption
  InitialStellarMass=Gal[galID].BulgeMass+Gal[galID].DiskMass;
  InitialBulgeMass=Gal[galID].BulgeMass;
  InitialDiskMass=Gal[galID].DiskMass;


  CentralVelDisp=Gal[galID].CentralVelDisp;
  CentralMass=2.0*CentralVelDisp*CentralVelDisp*Gal[galID].MergRadius/(G);

  //compute disruption radius
  DisruptionRadius=sqrt((0.5*Gal[galID].Vvir*Gal[galID].Vvir*Gal[galID].MergRadius*Gal[galID].MergRadius) /
                   (((G*CentralMass)/(2.0*Gal[galID].MergRadius))+(CentralVelDisp*CentralVelDisp)));

  if(DisruptionRadius < 1e-6)
 	  DisruptionRadius = 0.0;

  //disk disruption
  if(InitialDiskMass > 0.0 && DisruptionRadius < Gal[galID].StellarDiskRadius)
  {
	  Rd=(Gal[galID].StellarDiskRadius / 3.0);
	  Gal[galID].DiskMass = InitialDiskMass*(1.0-((1.0+DisruptionRadius/Rd)*exp(-DisruptionRadius/Rd)));
	  Gal[galID].StellarDiskRadius = DisruptionRadius;

	  if(Gal[galID].DiskMass > InitialDiskMass)
	  	Gal[galID].DiskMass=InitialDiskMass;

	  if(Gal[galID].DiskMass < 1e-6)
	  	Gal[galID].DiskMass=Gal[galID].StellarDiskRadius=0.0;

	  Gal[galID].MetalsDiskMass*=Gal[galID].DiskMass/InitialStellarMass;
  }


//bulge disruption
  //Re=pow(10.0,((log10(Gal[galID].Vvir/sqrt(2.0))-2.58)/0.21))/1000.0;
  if(Gal[galID].BulgeMass > 0.0 && DisruptionRadius < Gal[galID].BulgeSize)
  {
	  BulgeVelDisp=sqrt(G*(InitialBulgeMass)/(2.0*Gal[galID].BulgeSize));
	  Re=pow(10.0,((log10(BulgeVelDisp)-2.58)/0.21))/1000.0;
	  a=0.56*Re;
	  Gal[galID].BulgeMass=InitialBulgeMass *((DisruptionRadius/a)*(DisruptionRadius/a))
                                           /((1+(DisruptionRadius/a))*(1+(DisruptionRadius/a)));
	  Gal[galID].BulgeSize=DisruptionRadius;

	  if(Gal[galID].BulgeMass > InitialBulgeMass)
		  Gal[galID].BulgeMass = InitialBulgeMass;

	  if(Gal[galID].BulgeMass < 1e-6)
		  Gal[galID].BulgeMass=Gal[galID].BulgeSize=0.0;

	  Gal[galID].MetalsBulgeMass*=Gal[galID].BulgeMass/InitialBulgeMass;
  }

  if(Gal[galID].BulgeMass+Gal[galID].DiskMass < InitialStellarMass)
  {
	  Gal[Gal[galID].CentralGal].ICM+=(InitialStellarMass-Gal[galID].DiskMass-Gal[galID].BulgeMass);
#ifndef  POST_PROCESS_MAGS
	  sub_to_luminosities(galID, (Gal[galID].DiskMass+Gal[galID].BulgeMass)/InitialStellarMass*1.0);
#endif
  }

   //COLD GAS DISRUPTION
   //take cold gas mass into the central galaxy in the same proportion as diskmass
  if((InitialDiskMass-Gal[galID].DiskMass) > 0.0)
    {
	  LostFract=1.-Gal[galID].DiskMass/InitialDiskMass;
	  Gal[Gal[galID].CentralGal].HotGas+=Gal[galID].ColdGas*LostFract;
	  Gal[Gal[galID].CentralGal].MetalsHotGas+=Gal[galID].MetalsColdGas*LostFract;
	  Gal[galID].ColdGas*=(1.-LostFract);
	  Gal[galID].MetalsColdGas*=(1.-LostFract);
    }

  if (DisruptionRadius<Gal[galID].GasDiskRadius)
	  Gal[galID].GasDiskRadius=DisruptionRadius;


}



float erff(float x)
{
	float gammp(float a, float x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}


//If MCMC is defined this functions are already defined elsewhere
#ifndef MCMC

float gammp(float a, float x)
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



float gammq(float a, float x)
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



float gser(float a, float x)
{
  float gln;
  float sum, del, ap;
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




float gcf(float a, float x)
{
  float gln;
  int i;
  float an, b, c, d, del, h;

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

float gammpapprox(float a, float x, int psig)
{
  float y[18] = { 0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
    0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
    0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
    0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
    0.87126389619061517, 0.95698180152629142
  };
  float w[18] = { 0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
    0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
    0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
    0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
    0.085346685739338721, 0.085983275670394821
  };
  int j, ngau = 18;
  float xu, t, sum, ans;
  float a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
  float gln;

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




float gammln(float xx)
{
  int j;
  float x, tmp, y = 0, ser;

  float cof[14] = { 57.1562356658629235, -59.5979603554754912,
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



