#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

#ifdef UPDATETYPETWO

void update_type_two_coordinate_and_velocity(int tree, int i, int centralgal)
{
	int j, p;
	float tmppos;
	double Scale_V, dv;
	p=i;


#ifdef GUO10
    if(HaloGal[i].Type == 2)  /* Update positions of type 2's */
#else
    if(Gal[i].Type == 2)  /* Update positions of type 2's */
#endif
      {
#ifdef GUO10
    	int snapnum = HaloGal[i].SnapNum;
#else
    	int snapnum = Gal[i].SnapNum;
#endif
        Nids = CountIDs_snaptree[snapnum * Ntrees + tree];
        OffsetIDs = OffsetIDs_snaptree[snapnum * Ntrees + tree];

        size_t header_offset = 4 * sizeof(int) + 2 * TotSnaps * sizeof(int) + 2 * TotSnaps * Ntrees * sizeof(int) + 2
            * sizeof(int) * NtotHalos;

        IdList = (long long *) (TreeAuxData + header_offset);
        PosList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long));
        VelList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long) + TotIds * 3 * sizeof(float));

        IdList += OffsetIDs;
        PosList += 3 * OffsetIDs;
        VelList += 3 * OffsetIDs;
#ifdef GUO10
        get_coordinates(HaloGal[i].Pos, HaloGal[i].Vel, HaloGal[i].MostBoundID, tree, HaloGal[i].HaloNr,
        		HaloGal[i].SnapNum);
#else
        get_coordinates(Gal[i].Pos, Gal[i].Vel, Gal[i].MostBoundID, tree, Gal[i].HaloNr,
                		Gal[i].SnapNum);
#endif

//#ifdef SCALE_COSMOLOGY
#ifndef GUO10
        for(j = 0; j < 3; j++)
        	Gal[i].Pos[j] *= ScalePos;
#endif
//#endif

#ifdef GUO10
        for(j = 0; j < 3; j++)
          {
        	tmppos = wrap(-HaloGal[p].MergCentralPos[j] + HaloGal[p].Pos[j],BoxSize);
        	tmppos *=  sqrt(HaloGal[p].MergTime/HaloGal[p].OriMergTime);

        	HaloGal[p].Pos[j]=HaloGal[p].MergCentralPos[j] + tmppos;

        	if(HaloGal[p].Pos[j] < 0)
        		HaloGal[p].Pos[j] = BoxSize + HaloGal[p].Pos[j];
        	if(HaloGal[p].Pos[j] > BoxSize)
        		HaloGal[p].Pos[j] = HaloGal[p].Pos[j] - BoxSize;
          }
#else
        for(j = 0; j < 3; j++)
        {
        	tmppos = wrap(-Gal[p].MergCentralPos[j] + Gal[p].Pos[j],BoxSize);
#ifdef GUO13
        	tmppos *=  sqrt(Gal[p].MergTime/Gal[p].OriMergTime);
#else
        	tmppos *=  (Gal[p].MergTime/Gal[p].OriMergTime);
#endif
        	Gal[p].Pos[j]=Gal[p].MergCentralPos[j] + tmppos;

        	if(Gal[p].Pos[j] < 0)
        		Gal[p].Pos[j] = BoxSize + Gal[p].Pos[j];
        	if(Gal[p].Pos[j] > BoxSize)
        		Gal[p].Pos[j] = Gal[p].Pos[j] - BoxSize;
        }
#endif


#ifndef GUO10
//#ifdef SCALE_COSMOLOGY
	//add by Qi. 06/04/2012 to account for the scale of velocity field
        Scale_V = scale_v_cen(Halo[Gal[centralgal].HaloNr].SnapNum);
	 for (j = 0; j < 3 ; j++)
	{
		//printf("vel0=%0.3f cvel0=%0.3f ",Gal[p].Vel[j], Gal[centralgal].Vel[j]);
		dv = Gal[p].Vel[j] - Gal[centralgal].Vel[j]/Scale_V;
		dv *=sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[Gal[centralgal].HaloNr].SnapNum]/AA[Halo[Gal[centralgal].HaloNr].SnapNum]);
		Gal[p].Vel[j] = Gal[centralgal].Vel[j] + dv;
		//printf("vel1=%0.3f cel1=%0.3f ",Gal[p].Vel[j], Gal[centralgal].Vel[j]);
	}
	//printf("\n");
//#endif
#endif
      }
}

void get_coordinates(float *pos, float *vel, long long ID, int tree, int halonr, int snapnum)
{
  int m, k, start, nids;

  start = OffsetIDs_halo[TreeFirstHalo[tree] + halonr] - OffsetIDs;
  nids = CountIDs_halo[TreeFirstHalo[tree] + halonr];

  while(nids > 0)
  {
  	m = nids / 2;
  	if(IdList[start + m] == ID)
  	{
  		for(k = 0; k < 3; k++)
  		{
  			pos[k] = PosList[3 * (start + m) + k];
	      vel[k] = sqrt(AA[snapnum]) * VelList[3 * (start + m) + k];	/* to convert to peculiar velocity */
	    }

  		if(pos[0] == 0 && pos[1] == 0 && pos[2] == 0)
  		{
  			terminate	("This treeaux-files does not (yet) contain the coordinates\n for the desired output time!\n");
  		}

  		return;
  	}

  	if(IdList[start + m] < ID)
  	{
  		nids -= m;
  		start += m;
  	}
  	else
  	{
  		nids = m;
  	}
  }

  terminate("ID not found! - What's going on?");
  return;
}


/**@brief If USE_MEMORY_TO_MINIMIZE_IO ON - Routine to read in all
 *        the tree_aux data for one file into a pointer (ptr_auxdata) -
 *        later myfread() will pass the data onto each tree structure*/

void load_all_auxdata(int filenr)
{
  char buf[1000];
  FILE *fd;
  struct stat filestatus;


#ifndef MRII
  sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastDarkMatterSnapShot, filenr);
#else
  sprintf(buf, "%s/treedata/treeaux_sf1_%03d.%d", SimulationDir, LastDarkMatterSnapShot, filenr);
#endif


  if(stat(buf, &filestatus) != 0)                  /* seems not to exist */
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  size_t bytes = filestatus.st_size;

  TreeAuxData = mymalloc("TreeAuxData", bytes); 
  myfread(TreeAuxData, 1, bytes, fd);
  fclose(fd);
  
  int *header = TreeAuxData;

  NtotHalos = header[0];
  TotIds = header[1];
  Ntrees = header[2];
  TotSnaps = header[3];

  CountIDs_snaptree = header + 4 + 2 * TotSnaps;
  OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * Ntrees;
  CountIDs_halo = OffsetIDs_snaptree + TotSnaps * Ntrees;
  OffsetIDs_halo = CountIDs_halo + NtotHalos;
}




#endif

