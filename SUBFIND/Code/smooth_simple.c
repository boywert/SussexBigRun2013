#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"

#if (defined(DIVBCLEANING_DEDNER) || defined(SMOOTH_ROTB) || defined(BSMOOTH) || defined(SCAL_PRO_CLEAN)) || defined(VECT_POTENTIAL) || (defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS)

#if defined(LT_STELLAREVOLUTION)
int smooth_isactive(int);
#endif


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct smoothdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
#if (defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT))
/* #if defined(LT_SEvDbg) */
/*   MyIDType ID; */
/* #endif */
  int Type;
#endif
}
 *SmoothDataIn, *SmoothDataGet;


static struct smoothdata_out
{
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif
#ifdef VECT_POTENTIAL
  MyFloat SmoothA[3];
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
  MyFloat SmoothDivB;
#endif
#ifdef SMOOTH_ROTB
  MyFloat SmoothRotB[3];
#endif
#if defined(BSMOOTH)
  MyFloat BSmooth[3];
#endif
#if defined(BLACK_HOLES)
  MyLongDouble SmoothedEntr;
#endif

#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
  MyFloat Zsmooth[LT_NMetP];
#else
  float Zsmooth;
#endif
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
  float SmoothDens;
  int SmoothNgb;
#endif
#endif
#if defined(LT_SMOOTH_XCLD)
  float XCLDsmooth;
#endif
#if defined(LT_TRACK_WINDS)
  float AvgHsml;
#endif
  
  MyFloat DensityNorm;

}
 *SmoothDataResult, *SmoothDataOut;

void smoothed_values(void)
{
  int ngrp, sendTask, recvTask, place, nexport, nimport;
  int i, j, ndone, ndone_flag, dummy;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;

#if defined(BSMOOTH)
  int Smooth_Flag = 0;
  double dB[3];
#endif

#if defined(LT_SMOOTH_Z) && defined(LT_SMOOTH_ALLMETALS)
  int sp;
#endif
#ifdef LT_SMOOTH_SIZE
  double SmoothSize;
#endif

  /* Display information message that this step is executed on Task 0 ... */
  if(ThisTask == 0)
    {
      printf("Updating SPH interpolants for:"
#if defined(SMOOTH_PHI)
	     " (Phi - Dedner)"
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
	     " (Vect A)"
#endif /* VECT_POTENTIAL */
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
	     " (DivB) "
#endif
#ifdef SMOOTH_ROTB
	     " (rotB)"
#endif /* SMOOTH_ROTB */
#if defined(BSMOOTH)
	     " (B)"
#endif /* BSMOOTH */
#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
	     " (rho around stars)"
#endif
#ifdef LT_SMOOTH_Z
	     " (metallicity)"
#endif
#ifdef LT_SMOOTH_XCLD
	     " (cloud fraction)"
#endif
#ifdef LT_TRACK_WINDS
	     " (Hsml)"
#endif
	     "\n");
#ifdef BSMOOTH
      printf("Flag_FullStep = %d, Main TimestepCounts = %d\n", Flag_FullStep, All.MainTimestepCounts);
#endif
    }
#if defined(BSMOOTH)
  if(Flag_FullStep == 1)
    {
      if((All.MainTimestepCounts % All.BSmoothInt == 0) && (All.BSmoothInt >= 0))
	{
	  Smooth_Flag = 1;
	  if(ThisTask == 0)
	    printf("Smoothing B %d, %f\n", All.BSmoothInt, All.BSmoothFrac);
	}
      All.MainTimestepCounts++;
    }
#endif

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct smoothdata_in) + sizeof(struct smoothdata_out) +
					     sizemax(sizeof(struct smoothdata_in),
						     sizeof(struct smoothdata_out))));

  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
  unsigned int appended, tot_appended, alreadyactive, tot_alreadyactive;

  appended = append_chemicallyactive_particles(&alreadyactive);
  MPI_Reduce(&appended, &tot_appended, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&alreadyactive, &tot_alreadyactive, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0 && (tot_appended > 0 || tot_alreadyactive > 0))
    printf("%u chemically active particles queued for smoothing calculation (%u already active)..\n",
	   tot_appended, tot_alreadyactive);
  fflush(stdout);
#endif

  CPU_Step[CPU_SMTHMISC] += measure_time();
  t0 = second();

  i = FirstActiveParticle;	/* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	{
#if !defined(LT_STELLAREVOLUTION)
	  if(density_isactive(i))
#else
	  if(smooth_isactive(i))
#endif
	    {
	      if(smoothed_evaluate(i, 0, &nexport, Send_count) < 0)
		break;
	    }
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      tstart = second();

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timewait1 += timediff(tstart, tend);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      SmoothDataGet =
	(struct smoothdata_in *) mymalloc("SmoothDataGet", nimport * sizeof(struct smoothdata_in));
      SmoothDataIn =
	(struct smoothdata_in *) mymalloc("SmoothDataIn", nexport * sizeof(struct smoothdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  SmoothDataIn[j].Pos[0] = P[place].Pos[0];
	  SmoothDataIn[j].Pos[1] = P[place].Pos[1];
	  SmoothDataIn[j].Pos[2] = P[place].Pos[2];
	  SmoothDataIn[j].Hsml = PPP[place].Hsml;
#if defined(LT_STELLAREVOLUTION)
/* #if defined(LT_SEvDbg) */
/* 	  SmoothDataIn[j].ID = P[place].ID; */
/* #endif */
#if !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
	  SmoothDataIn[j].Type = P[place].Type;
#endif
#endif
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
#ifdef LT_SMOOTH_NGB
          SmoothDataIn[j].SmoothHsml = SphP[place].SmoothHsml;
#else
          if((SmoothDataIn[j].SmoothHsml =
              SmoothDataIn[j].Hsml * All.SmoothRegionSize) > All.SmoothRegionSizeMax)
            SmoothDataIn[j].SmoothHsml = All.SmoothRegionSizeMax;
#endif
#endif
	  memcpy(SmoothDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&SmoothDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct smoothdata_in), MPI_BYTE,
			       recvTask, TAG_SMOOTH_A,
			       &SmoothDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct smoothdata_in), MPI_BYTE,
			       recvTask, TAG_SMOOTH_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(SmoothDataIn);
      SmoothDataResult =
	(struct smoothdata_out *) mymalloc("SmoothDataResult", nimport * sizeof(struct smoothdata_out));
      SmoothDataOut =
	(struct smoothdata_out *) mymalloc("SmoothDataOut", nexport * sizeof(struct smoothdata_out));


      /* now do the particles that were sent to us */

      tstart = second();
      for(j = 0; j < nimport; j++)
	smoothed_evaluate(j, 1, &dummy, &dummy);
      tend = second();
      timecomp2 += timediff(tstart, tend);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);


      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&SmoothDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct smoothdata_out),
			       MPI_BYTE, recvTask, TAG_SMOOTH_B,
			       &SmoothDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct smoothdata_out),
			       MPI_BYTE, recvTask, TAG_SMOOTH_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm2 += timediff(tstart, tend);


      /* add the result to the local particles */
      tstart = second();
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  if(P[place].Type == 0)
	    {
#if defined(SMOOTH_PHI)
	      SphP[place].SmoothPhi += SmoothDataOut[j].SmoothPhi;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
	      SphP[place].SmoothA[0] += SmoothDataOut[j].SmoothA[0];
	      SphP[place].SmoothA[1] += SmoothDataOut[j].SmoothA[1];
	      SphP[place].SmoothA[2] += SmoothDataOut[j].SmoothA[2];
#endif /* VECT_POTENTIAL */
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
	      SphP[place].SmoothDivB += SmoothDataOut[j].SmoothDivB;
#endif
#ifdef SMOOTH_ROTB
	      SphP[place].SmoothedRotB[0] += SmoothDataOut[j].SmoothRotB[0];
	      SphP[place].SmoothedRotB[1] += SmoothDataOut[j].SmoothRotB[1];
	      SphP[place].SmoothedRotB[2] += SmoothDataOut[j].SmoothRotB[2];
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
	      SphP[place].BSmooth[0] += SmoothDataOut[j].BSmooth[0];
	      SphP[place].BSmooth[1] += SmoothDataOut[j].BSmooth[1];
	      SphP[place].BSmooth[2] += SmoothDataOut[j].BSmooth[2];
#endif /* BSMOOTH */
#ifdef LT_SMOOTH_Z
#ifdef LT_SMOOTH_ALLMETALS
              for(sp = 0; sp < LT_NMetP; sp++)
                SphP[place].Zsmooth[sp] += SmoothDataOut[j].Zsmooth[sp];
#else
	      SphP[place].Zsmooth += SmoothDataOut[j].Zsmooth;
#endif
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
	      SphP[place].SmoothDens += SmoothDataOut[j].SmoothDens;
	      SphP[place].SmoothNgb += SmoothDataOut[j].SmoothNgb;
#endif
#endif /* LT_SMOOTH_Z */

#ifdef LT_SMOOTH_CLDX
	      SphP[place].XCLDsmooth += SmoothDataOut[j].XCLDsmooth;
#endif /* LT_SMOOTH_CLDX */
#ifdef LT_TRACK_WINDS
	      SphP[place].AvgHsml += SmoothDataOut[j].AvgHsml;
#endif /* LT_TRACK_WINDS */
	      SphP[place].DensityNorm += SmoothDataOut[j].DensityNorm;
	    }

#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
	  if(P[place].Type == 4)
	    MetP[P[place].pt.MetID].weight += SmoothDataOut[j].DensityNorm;
#endif
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);


      myfree(SmoothDataOut);
      myfree(SmoothDataResult);
      myfree(SmoothDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);



  /* do final operations on results */
  tstart = second();
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#if !defined(LT_STELLAREVOLUTION)
      if(density_isactive(i))
#else
      if(smooth_isactive(i))
#endif
        {
          if(P[i].Type == 0)
            {
#if defined(SMOOTH_PHI)
              SphP[i].SmoothPhi /= SphP[i].DensityNorm;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
              SphP[i].SmoothA[0] /= SphP[i].DensityNorm;
              SphP[i].SmoothA[1] /= SphP[i].DensityNorm;
              SphP[i].SmoothA[2] /= SphP[i].DensityNorm;
              SphP[i].APred[0] = SphP[i].SmoothA[0];
              SphP[i].APred[1] = SphP[i].SmoothA[1];
              SphP[i].APred[2] = SphP[i].SmoothA[2];
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
              SphP[i].SmoothDivB /= SphP[i].DensityNorm;
#endif
#ifdef SMOOTH_ROTB
              SphP[i].SmoothedRotB[0] /= SphP[i].DensityNorm;
              SphP[i].SmoothedRotB[1] /= SphP[i].DensityNorm;
              SphP[i].SmoothedRotB[2] /= SphP[i].DensityNorm;
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
              SphP[i].BSmooth[0] /= SphP[i].DensityNorm;
              SphP[i].BSmooth[1] /= SphP[i].DensityNorm;
              SphP[i].BSmooth[2] /= SphP[i].DensityNorm;

              if(Smooth_Flag == 1)
                {
                  dB[0] = All.BSmoothFrac * (SphP[i].BSmooth[0] - SphP[i].b2.BPred[0]);
                  dB[1] = All.BSmoothFrac * (SphP[i].BSmooth[1] - SphP[i].b2.BPred[1]);
                  dB[2] = All.BSmoothFrac * (SphP[i].BSmooth[2] - SphP[i].b2.BPred[2]);
                  SphP[i].b2.BPred[0] += dB[0];
                  SphP[i].b2.BPred[1] += dB[1];
                  SphP[i].b2.BPred[2] += dB[2];
#ifndef EULERPOTENTIALS
                  SphP[i].b1.B[0] += dB[0];
                  SphP[i].b1.B[1] += dB[1];
                  SphP[i].b1.B[2] += dB[2];
#endif
                }
#endif /* BSMOOTH */

#ifdef LT_SMOOTH_Z		/* > ============================================== < */
              /* >  SMOOTH Z                                      < */

#if !defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)
#if defined(LT_SMOOTH_ALLMETALS)
              for(sp = 0; sp < LT_NMetP; sp++)
                SphP[i].Zsmooth[sp] /= SphP[i].DensityNorm;
#else
              SphP[i].Zsmooth /= SphP[i].DensityNorm;
#endif
#else
              SphP[i].Zsmooth /= SphP[i].SmoothDens;
#endif
#endif /* LT_SMOOTH:Z */
#ifdef LT_SMOOTH_XCLD
              SphP[i].XCLDsmooth /= SphP[i].DensityNorm;
#endif
#ifdef LT_TRACK_WINDS
              SphP[i].AvgHsml /= SphP[i].DensityNorm;
#endif
#if defined(LT_SMOOTH_SIZE)
              AvgSmoothN++;

              if((SmoothSize = PPP[i].Hsml * All.SmoothRegionSize) > All.SmoothRegionSizeMax)
                SmoothSize = All.SmoothRegionSizeMax;
              AvgSmoothSize += SmoothSize;
              if(SmoothSize < MinSmoothSize)
                MinSmoothSize = SmoothSize;
              if(SmoothSize > MaxSmoothSize)
                MaxSmoothSize = SmoothSize;

              AvgSmoothNgb += SphP[i].SmoothNgb;
              if(SphP[i].SmoothNgb < MinSmoothNgb)
                MinSmoothNgb = SphP[i].SmoothNgb;
              if(SphP[i].SmoothNgb > MaxSmoothNgb)
                MaxSmoothNgb = SphP[i].SmoothNgb;
#endif
            }

        }
    }
  tend = second();
  timecomp1 += timediff(tstart, tend);

#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
  drop_chemicallyactive_particles();
#endif

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_SMTHCOMPUTE] += timecomp;
  CPU_Step[CPU_SMTHWAIT] += timewait;
  CPU_Step[CPU_SMTHCOMM] += timecomm;
  CPU_Step[CPU_SMTHMISC] += timeall - (timecomp + timewait + timecomm);

}



/*! This function represents the core of the SPH density computation. The
*  target particle may either be local, or reside in the communication
*  buffer.
*/
int smoothed_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, listindex = 0;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3, hinv4;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat *pos;
  double DensityNorm = 0;

#ifdef SMOOTH_PHI
  double SmoothPhi = 0.0;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
  double smoothA[3];

  smoothA[0] = smoothA[1] = smoothA[2] = 0.0;
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
  double SmoothDivB = 0.0;
#endif
#ifdef SMOOTH_ROTB
  double smoothrotb[3];
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
  double BSmooth[3];

  BSmooth[0] = BSmooth[1] = BSmooth[2] = 0;
#endif /* BSMOOTH */

/*  MyIDType myID; */

#if !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
  int Type;
#endif

#if defined(LT_SMOOTH_Z)	/* ======= LT_SMOOTH_Z */

#if defined(LT_SMOOTH_ALLMETALS)
  int sp;
  MyDouble Zrho[LT_NMetP];
#else
  MyDouble getmetallicity;
  MyDouble Zrho = 0;
#endif  
  
  double SmoothSize;

#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)	/* >>>>> */
  int SmoothNgb;
  double su;
  DOUBLE SmoothDens = 0;
  DOUBLE SmoothHsml, shinv, shinv3, shinv4;

#ifdef LT_SMOOTH_NGB
  int smoothcc;
#endif
#endif /* <<<<< closes #if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB) */

#endif /* closes LT_SMOOTH_Z */

#ifdef LT_SMOOTH_XCLD
  DOUBLE XCLDsmooth;
#endif

#ifdef LT_TRACK_WINDS
  DOUBLE AvgHsml;
#endif

#if defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
  double swk, sdwk;
#endif

#ifdef SMOOTH_ROTB
  smoothrotb[0] = smoothrotb[1] = smoothrotb[2] = 0;
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
  BSmooth[0] = BSmooth[1] = BSmooth[2] = 0;
#endif /* BSMOOTH */

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
#if defined(LT_STELLAREVOLUTION)
/* #if defined(LT_SEvDbg) */
/*       myID = P[target].ID; */
/* #endif */
#if !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
      Type = P[target].Type;
#endif
#endif
    }
  else
    {
      pos = SmoothDataGet[target].Pos;
      h = SmoothDataGet[target].Hsml;
#if defined(LT_STELLAREVOLUTION)
/* #if defined(LT_SEvDbg) */
/*       myID = SmoothDataGet[target].ID; */
/* #endif */
#if !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
      Type = SmoothDataGet[target].Type;
#endif
#endif
    }

  h2 = h * h;

  kernel_hinv(h, &hinv, &hinv3, &hinv4);

#if defined(LT_SMOOTH_Z) && defined(LT_SMOOTH_ALLMETALS)
  memset(Zrho, 0, LT_NMetP * sizeof(MyDouble));
#endif
  
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
#ifdef LT_SMOOTH_SIZE
  if((SmoothHsml = h * All.SmoothRegionSize) > All.SmoothRegionSizeMax)
    SmoothHsml = All.SmoothRegionSizeMax;
#endif
#ifdef LT_SMOOTH_NGB
  SmoothHsml = h;
  memset(NearestSmooth, 0, All.DesNumNbSmooth * sizeof(float));
#endif
  kernel_hinv(SmoothHsml, &shinv, &shinv3, &shinv4);
  SmoothNgb = 0;
#endif


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = SmoothDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifndef KD_FRICTION
          numngb_inbox = ngb_treefind_variable(&pos[0], h, target, &startnode, mode, nexport, nsend_local);
#else
          numngb_inbox = ngb_treefind_variable(&pos[0], h, target, &startnode, mode, nexport, nsend_local, 0);
#endif
	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < h2)
		{
		  r = sqrt(r2);

#ifdef LT_SMOOTH_NGB
		  for(smoothcc = 0; smoothcc < All.DesNumNgbSmooth; smoothcc++)
		    if(Nearest[smoothcc] > (float) r)
		      Nearest[smoothcc] = (float) r;
#endif
		  u = r * hinv;
		  kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
		  mass_j = P[j].Mass;
		  wk /= SphP[j].d.Density;

#if defined(SMOOTH_PHI)
		  SmoothPhi += mass_j * wk * SphP[j].PhiPred;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
		  smoothA[0] += mass_j * wk * SphP[j].APred[0];
		  smoothA[1] += mass_j * wk * SphP[j].APred[1];
		  smoothA[2] += mass_j * wk * SphP[j].APred[2];
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
		  SmoothDivB += mass_j * wk * SphP[j].divB;
#endif
#ifdef SMOOTH_ROTB
		  smoothrotb[0] += mass_j * wk * SphP[j].RotB[0];
		  smoothrotb[1] += mass_j * wk * SphP[j].RotB[1];
		  smoothrotb[2] += mass_j * wk * SphP[j].RotB[2];
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
		  BSmooth[0] += mass_j * wk * SphP[j].b2.BPred[0];
		  BSmooth[1] += mass_j * wk * SphP[j].b2.BPred[1];
		  BSmooth[2] += mass_j * wk * SphP[j].b2.BPred[2];
#endif /* BSMOOTH */

#if defined(LT_SMOOTH_Z)	/* ========== LT_SMOOTH_Z */

#if defined(LT_SMOOTH_ALLMETALS)
                  for(sp = 0; sp < LT_NMetP; sp++)
                    Zrho[sp] += (MyDouble) SphP[j].Metals[sp] / P[j].Mass * mass_j * wk;            
#else                  
#if !defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)	/* SMOOTH_SIZE and SMOOT_NGB are NOT defined */
		  getmetallicity = get_metallicity(j, -1);
		  Zrho += (MyDouble) getmetallicity * mass_j * wk;
                  swk = wk;
#else /* SMOOTH_SIZE or SMOOTH_NGB are defined */
		  if(r <= SmoothHsml)
		    {
		      SmoothNgb++;

		      su = r * shinv;

		      kernel_main(su, shinv3, shinv4, &swk, &sdwk, -1);

                      swk /= SphP[j].d.Density;
                      
		      getmetallicity = get_metallicity(j, -1);
		      Zrho += (MyDouble) getmetallicity * mass_j * swk;

		      SmoothDens += mass_j * swk;
                    }
#endif /* closes SMOOTH_SIZE || SMOOTH_NGB */
#endif /* closes LT_SMOOTH_ALLMETALS */                                                     
#endif /* ========== closes LT_SMOOTH_Z */

#ifdef LT_SMOOTH_XCLD
		  XCLDsmooth += (MyDouble) SphP[j].x * mass_j * wk;
#endif /* LT_SMOOTH_XCLD */

#ifdef LT_TRACK_WINDS
		  AvgHsml += (MyDouble) SphP[j].Hsml * mass_j * wk;
#endif /* LT_TRACK_WINDS */

		  DensityNorm += mass_j * wk;
                }
	    }
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = SmoothDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  if(mode == 0)
    {
#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
      if(Type == 4)
	MetP[P[target].pt.MetID].weight = DensityNorm;

      if(Type == 0)   /* protect to write into SphP with traget not to be a gas particle ! */
	{
#endif

#if defined(SMOOTH_PHI)
          SphP[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
          SphP[target].SmoothA[0] = smoothA[0];
          SphP[target].SmoothA[1] = smoothA[1];
          SphP[target].SmoothA[2] = smoothA[2];
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
          SphP[target].SmoothDivB = SmoothDivB;
#endif
#ifdef SMOOTH_ROTB
          SphP[target].SmoothedRotB[0] = smoothrotb[0];
          SphP[target].SmoothedRotB[1] = smoothrotb[1];
          SphP[target].SmoothedRotB[2] = smoothrotb[2];
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
          SphP[target].BSmooth[0] = BSmooth[0];
          SphP[target].BSmooth[1] = BSmooth[1];
          SphP[target].BSmooth[2] = BSmooth[2];
#endif /* BSMOOTH */

#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
          for(sp = 0; sp < LT_NMetP; sp++)
            SphP[target].Zsmooth[sp] = (float)Zrho[sp];
#else
          SphP[target].Zsmooth = (float) Zrho;
#endif
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
          SphP[target].SmoothDens = (float) SmoothDens;
          SphP[target].SmoothNgb = SmoothNgb;
#endif
#endif /* LT_SMOOTH_Z */

#if defined(LT_SMOOTH_XCLD)
          SphP[target].XCLDsmooth = (float) XCLDsmooth;
#endif /* LT_SMOOTH_XCLD */

#if defined(LT_TRACK_WINDS)
          SphP[target].AvgHsml = (float) AvgHsml;
#endif /* LT_TRACK_WINDS */

          SphP[target].DensityNorm = DensityNorm;

#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
        }
#endif

    }
  else
    {
#if defined(SMOOTH_PHI)
      SmoothDataResult[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */
#if defined(VECT_POTENTIAL)
      SmoothDataResult[target].SmoothA[0] = smoothA[0];
      SmoothDataResult[target].SmoothA[1] = smoothA[1];
      SmoothDataResult[target].SmoothA[2] = smoothA[2];
#endif
#if defined(DIVBCLEANING_DEDNER) || defined(SCAL_PRO_CLEAN)
      SmoothDataResult[target].SmoothDivB = SmoothDivB;
#endif
#ifdef SMOOTH_ROTB
      SmoothDataResult[target].SmoothRotB[0] = smoothrotb[0];
      SmoothDataResult[target].SmoothRotB[1] = smoothrotb[1];
      SmoothDataResult[target].SmoothRotB[2] = smoothrotb[2];
#endif /* SMOOTH_ROTB */

#if defined(BSMOOTH)
      SmoothDataResult[target].BSmooth[0] = BSmooth[0];
      SmoothDataResult[target].BSmooth[1] = BSmooth[1];
      SmoothDataResult[target].BSmooth[2] = BSmooth[2];
#endif /* BSMOOTH */

#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
      for(sp = 0; sp < LT_NMet; sp++)
        SmoothDataResult[target].Zsmooth[sp] = (float)Zrho[sp];
#else
      SmoothDataResult[target].Zsmooth = (float) Zrho;
#endif      
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
      SmoothDataResult[target].SmoothDens = (float) SmoothDens;
      SmoothDataResult[target].SmoothNgb = SmoothNgb;
#endif
#endif /* LT_SMOOTH_Z */

#if defined(LT_SMOOTH_XCLD)
      SmoothDataResult[target].XCLDsmooth = (float) XCLDsmooth;
#endif /* LT_SMOOTH_XCLD */

#if defined(LT_TRACK_WINDS)
      SmoothDataResult[target].AvgHsml = (float) AvgHsml;
#endif /* LT_TRACK_WINDS */
      SmoothDataResult[target].DensityNorm = DensityNorm;
      
    }

  return 0;
}


#if defined(LT_STELLAREVOLUTION)

int smooth_isactive(int i)
{
  if(P[i].TimeBin < 0)
    return 0;

#if (defined(DIVBCLEANING_DEDNER) || defined(SMOOTH_ROTB) || defined(BSMOOTH) || defined(SCAL_PRO_CLEAN)) || defined(VECT_POTENTIAL) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS)
  if(P[i].Type == 0)
    return 1;
#endif

#if defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)
  if((P[i].Type & 15) == 4)
    if(TimeBinActive[MetP[P[i].pt.MetID].ChemTimeBin])
      return 1;
#endif

#ifdef BLACK_HOLES
  if(P[i].Type == 5)
    return 1;
#endif

  return 0;

}
#endif


#endif
