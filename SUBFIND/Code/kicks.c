#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

void apply_long_range_kick(int, int);

void do_first_halfstep_kick(void)
{
  int i;
  int ti_step, tstart, tend;

#ifdef PMGRID
  if(All.PM_Ti_begstep == All.Ti_Current)       /* need to do long-range kick */
    {
      ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      tstart = All.PM_Ti_begstep;
      tend = tstart + ti_step / 2;

      apply_long_range_kick(tstart, tend);
    }
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
     {
       ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;

       tstart = P[i].Ti_begstep; /* beginning of step */
       tend = P[i].Ti_begstep + ti_step / 2;     /* midpoint of step */

       do_the_kick(i, tstart, tend, P[i].Ti_current);
     }
}

void do_second_halfstep_kick(void)
{
  int i, j;
  int ti_step, tstart, tend;

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)       /* need to do long-range kick */
    {
      ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      tstart = All.PM_Ti_begstep + ti_step / 2;
      tend = tstart + ti_step / 2;

      apply_long_range_kick(tstart, tend);
    }
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
     {
       ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;

       tstart = P[i].Ti_begstep + ti_step / 2;    /* midpoint of step */
       tend = P[i].Ti_begstep + ti_step; /* end of step */

       do_the_kick(i, tstart, tend, P[i].Ti_current);

       /* after completion of a full step, set the predication values of SPH quantities
        * to the current values. They will then predicted along in drift operations
        */
       if(P[i].Type == 0)
         {
           for(j=0;j<3 ;j++)
             SphP[i].VelPred[j] = P[i].Vel[j];

           SphP[i].EntropyPred = SphP[i].Entropy;

#ifdef EOS_DEGENERATE
	   for(j = 0; j < 3; j++)
	     SphP[i].xnucPred[j] = SphP[i].xnuc[j];
#endif

           SphP[i].Pressure = get_pressure(i);

           set_predicted_sph_quantities_for_extra_physics(i);
         }
     }
}


#ifdef PMGRID
void apply_long_range_kick(int tstart, int tend)
{
  int i, j;
  double dt_gravkick, dvel[3];

  if(All.ComovingIntegrationOn)
    dt_gravkick = get_gravkick_factor(tstart, tend);
  else
    dt_gravkick = (tend - tstart) * All.Timebase_interval;

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)      /* do the kick, only collisionless particles */
        {
          dvel[j] = P[i].GravPM[j] * dt_gravkick;
          P[i].Vel[j] += dvel[j];
          P[i].dp[j] += P[i].Mass * dvel[j];
        }
#ifdef DISTORTIONTENSORPS
      do_distortion_tensor_kick(i, dt_gravkick);
#endif
    }
}
#endif


void do_the_kick(int i, int tstart, int tend, int tcurrent)
{
  int j;
  double dv[3];
  double dt_entr, dt_gravkick, dt_hydrokick;

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
     }
  else
    {
      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
     }


  /* do the kick */
  for(j = 0; j < 3; j++)
    {
      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;
#ifdef RELAXOBJECT
      dv[j] -= P[i].Vel[j] * All.RelaxFac * dt_gravkick;
#endif

      P[i].Vel[j] += dv[j];
      P[i].dp[j] += P[i].Mass * dv[j];
    }

#ifdef DISTORTIONTENSORPS
  do_distortion_tensor_kick(i, dt_gravkick);
#endif


  if(P[i].Type == 0)   /* kick for SPH quantities */
    {
      for(j = 0; j < 3; j++)
	{
	  dv[j] = SphP[i].a.HydroAccel[j] * dt_hydrokick;
	  P[i].Vel[j] += dv[j];
	  P[i].dp[j] += P[i].Mass * dv[j];
	}

      double dEntr = SphP[i].e.DtEntropy * dt_entr;

#if defined(EOS_DEGENERATE)
      dEntr *= All.UnitTime_in_s;
#endif

      SphP[i].Entropy = DMAX(SphP[i].Entropy + dEntr, 0.5 * SphP[i].Entropy);

      check_particle_for_temperature_minimum(i);

      do_sph_kick_for_extra_physics(i, tstart, tend, dt_entr);
    }
}

void set_predicted_sph_quantities_for_extra_physics(int i)
{
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS) && !defined(VECT_POTENTIAL)
  int j1;
  for(j1=0;j1<3 ;j1++)
    SphP[i].b2.BPred[j1] = SphP[i].b1.B[j1];
#endif
#ifdef VECT_POTENTIAL
  int j2;
  for(j2=0;j2<3 ;j2++)
    SphP[i].APred[j2] = SphP[i].A[j2];
#endif
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
   SphP[i].PhiPred = SphP[i].Phi;
#endif
}


void do_sph_kick_for_extra_physics(int i, int tstart, int tend, double dt_entr)
{
  int j;

#ifdef MAGNETIC
  double dt_mag; 
  if(All.ComovingIntegrationOn)
    dt_mag = get_magkick_factor(tstart, tend);
  else
    dt_mag = (tend - tstart) * All.Timebase_interval;
#endif

  for(j = 0; j < 3; j++)
    {
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS) && !defined(VECT_POTENTIAL)
      SphP[i].b1.B[j] += SphP[i].DtB[j] * dt_mag;
#endif
#ifdef VECT_POTENTIAL
      SphP[i].A[j] += SphP[i].DtA[j] * dt_entr;
#endif
    }

#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
  SphP[i].Phi += SphP[i].DtPhi * dt_mag;
#endif
#ifdef TIME_DEP_ART_VISC
  SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
  SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
  if(SphP[i].alpha < All.AlphaMin)
  SphP[i].alpha = All.AlphaMin;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
  SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
#ifdef VORONOI_RELAX_VIA_VISC
  if(SphP[i].alpha < All.ArtBulkViscConst / 128.0 / 128.0)
  SphP[i].alpha = All.ArtBulkViscConst / 128.0 / 128.0;
#else
  if(SphP[i].alpha < All.ArtBulkViscConst / 128.0)
  SphP[i].alpha = All.ArtBulkViscConst / 128.0;
#endif
#endif
#ifdef TIME_DEP_MAGN_DISP
  SphP[i].Balpha += SphP[i].DtBalpha * dt_entr;
  SphP[i].Balpha = DMIN(SphP[i].Balpha, All.ArtMagDispConst);
  if(SphP[i].Balpha < All.ArtMagDispMin)
    SphP[i].Balpha = All.ArtMagDispMin;
#endif

#ifdef NUCLEAR_NETWORK
  for(j = 0; j < EOS_NSPECIES; j++)
    SphP[i].xnuc[j] += SphP[i].dxnuc[j] * dt_entr * All.UnitTime_in_s;

  network_normalize(SphP[i].xnuc, &SphP[i].Entropy, &All.nd, &All.nw);
#endif

#ifdef CHEMISTRY
  /* update the chemical abundances for the new density and temperature */
  double a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
  double a_end = All.TimeBegin * exp(tend * All.Timebase_interval);
  int mode;
  /* time in cosmic expansion parameter */
  compute_abundances(mode = 1, i, a_start, a_end);
#endif
}



#ifdef DISTORTIONTENSORPS
void do_distortion_tensor_kick(int i, double dt_gravkick)
  {
   /* for the distortion 'velocity part', so only the lower two 3x3 submatrices will be != 0 */
  MyDouble dv_distortion_tensorps[6][6];
  int j, j1, j2;

     /* now we do the distortiontensor kick */
  for(j1 = 0; j1 < 3; j1++)
    for(j2 = 0; j2 < 3; j2++)
      {
        dv_distortion_tensorps[j1 + 3][j2] = 0.0;
        dv_distortion_tensorps[j1 + 3][j2 + 3] = 0.0;

        /* the 'acceleration' is given by the product of tidaltensor and distortiontensor */
        for(j = 0; j < 3; j++)
          {
            dv_distortion_tensorps[j1 + 3][j2] +=
              P[i].tidal_tensorps[j1][j] * P[i].distortion_tensorps[j][j2];
            dv_distortion_tensorps[j1 + 3][j2 + 3] +=
              P[i].tidal_tensorps[j1][j] * P[i].distortion_tensorps[j][j2 + 3];
          }
        dv_distortion_tensorps[j1 + 3][j2] *= dt_gravkick;
        dv_distortion_tensorps[j1 + 3][j2 + 3] *= dt_gravkick;

        /* add it to the distortiontensor 'velocities' */
        P[i].distortion_tensorps[j1 + 3][j2] += dv_distortion_tensorps[j1 + 3][j2];
        P[i].distortion_tensorps[j1 + 3][j2 + 3] += dv_distortion_tensorps[j1 + 3][j2 + 3];
      }
  }
#endif
