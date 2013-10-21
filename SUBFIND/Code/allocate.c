#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* This routine allocates memory for
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  size_t bytes;

  double bytes_tot = 0;

  int NTaskTimesThreads;

  NTaskTimesThreads = maxThreads * NTask;

  Exportflag = (int *) mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);

#ifdef VORONOI
  Mesh_Send_count = (int *) mymalloc("Mesh_Send_count", sizeof(int) * NTask);
  Mesh_Send_offset = (int *) mymalloc("Mesh_Send_offset", sizeof(int) * NTask);
  Mesh_Recv_count = (int *) mymalloc("Mesh_Recv_count", sizeof(int) * NTask);
  Mesh_Recv_offset = (int *) mymalloc("Mesh_Recv_offset", sizeof(int) * NTask);
#endif


  ProcessedFlag = (unsigned char *) mymalloc("ProcessedFlag", bytes = All.MaxPart * sizeof(unsigned char));
  bytes_tot += bytes;

  NextActiveParticle = (int *) mymalloc("NextActiveParticle", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  NextInTimeBin = (int *) mymalloc("NextInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  PrevInTimeBin = (int *) mymalloc("PrevInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;


  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc("P", bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }

  if(All.MaxPartSph > 0)
    {
      bytes_tot = 0;

      if(!
	 (SphP =
	  (struct sph_particle_data *) mymalloc("SphP", bytes =
						All.MaxPartSph * sizeof(struct sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }

#ifdef LT_STELLAREVOLUTION
  if(All.MaxPartMet > 0)
    {
      bytes_tot = 0;

      if(!
	 (MetP =
	  (struct met_particle_data *) mymalloc("MetP", bytes =
						All.MaxPartMet * sizeof(struct met_particle_data))))
	{
	  printf("failed to allocate memory for `MetP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of MET data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
#endif

#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
  if(All.MaxPartBH > 0)
    {
      bytes_tot = 0;

      if(!
	 (BHP =
	  (struct bh_particle_data *) mymalloc("BHP", bytes =
						All.MaxPartBH * sizeof(struct bh_particle_data))))
	{
	  printf("failed to allocate memory for `BHP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of BH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
#endif

#ifdef WRITE_KEY_FILES
  if(All.MaxPart > 0)
    {
      bytes_tot = 0;
      KeyIndex = mymalloc("KeyIndex", bytes = All.MaxPart * sizeof(peanokey));
      bytes_tot += bytes;
      NPartPerKey = mymalloc("KeyNpart", bytes = All.MaxPart * sizeof(int));
      bytes_tot += bytes;
      PartKeyOffset = mymalloc("KeyOffset", bytes = All.MaxPart * sizeof(int));
      bytes_tot += bytes;
      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of KEY data.\n\n", bytes_tot / (1024.0 * 1024.0));
    }
#endif
}
