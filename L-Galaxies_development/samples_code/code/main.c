/** @file main.c
 * @brief The file containing the main() function for L-Galaxies.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "proto.h"

//Ntrees=4612
//tothalos=31151583

/**@file main.c
 * @brief Controlling function of L-Galaxies plus Construct Galaxies,
 *        Join Galaxies of progenitors, Evolve Galaxies and check
 *        Makefile options.
 * */
int main(int argc, char **argv) {
  int treenr, i, j, k, ll;
  int index;
  FILE *fa, *fb, *fc;
  char buf[1000];
  //int NFOFs=767;
  char MR_SampleName[1000];
  int *Unique_TreeList_MR, *Unique_FileList_MR, Unique_TreeCount_MR;
  long long *FOFIDs;
  int FileNumber_MR=250;
  int MaxTrees=10000;
#ifdef MRII
  char MRII_SampleName[1000];
  int *Unique_TreeList_MRII, *Unique_FileList_MRII, Unique_TreeCount_MRII;
  int FileNumber_MRII=70;
#endif
  int Total_Unique_TreeCount;
  int Unique_totNHalos, Unique_totNHalos_MR,*Unique_TreeNHalos, *Unique_TreeFirstHalo;
  int MAX_UNIQUE_TREES=40000;
  int *Aux_CountIDs_snaptree, *Aux_OffsetIDs_snaptree, *Aux_CountIDs_halo, *Aux_OffsetIDs_halo;
  long long Aux_TotIds;
  int *Aux_CountIDs_snap, *Aux_OffsetIDs_snap;
#ifdef MR
  int Nsnaps=64;
#endif
#ifdef MRII
  int Nsnaps=68;
#endif
#ifdef  HALO_MODEL
  int sample_type=4; //marcel_M200=4, marcel=3, cut_optimal=2, optimal=1
#else
  int sample_type=2; //marcel_M200=4, marcel=3, cut_optimal=2, optimal=1
#endif
  int snap_low=20, snap_high_MR=63, snap_high_MRII=67;

#ifdef  HALO_MODEL
  sprintf(MR_SampleName, "optimal_marcel_M200_meansample_allz_nh");
#ifdef MRII
  sprintf(MRII_SampleName, "optimal_wx1_MRIIsample_allz_nh");
#endif
#else
  sprintf(MR_SampleName, "cut_optimalsample_allz_nh");
 #ifdef MRII
   sprintf(MRII_SampleName, "optimal_MRIIsample_allz_nh");
 #endif
#endif
  /*Reads the parameter file, given as an argument at run time.*/
  read_parameter_file(argv[1]);

  //MR
  Unique_TreeList_MR = malloc(sizeof(int) * MAX_UNIQUE_TREES);
  Unique_FileList_MR = malloc(sizeof(int) * MAX_UNIQUE_TREES);
  Unique_TreeCount_MR=0;

  //FIND LIST OF UNIQUE LIST OF TREES AND FILES
  printf("\nFinding Unique Trees for MR\n");
  Unique_TreeCount_MR=find_unique_tree_list(Unique_TreeList_MR, Unique_FileList_MR, &Unique_TreeCount_MR,
  		                                      MR_SampleName, FileNumber_MR, snap_low, snap_high_MR, MAX_UNIQUE_TREES);
  //Unique_TreeCount_MR=1.;
  Total_Unique_TreeCount=Unique_TreeCount_MR;
#ifdef MRII
  printf("Finding Unique Trees for MRII\n");
  Unique_TreeList_MRII = malloc(sizeof(int) * MAX_UNIQUE_TREES);
  Unique_FileList_MRII = malloc(sizeof(int) * MAX_UNIQUE_TREES);
  Unique_TreeCount_MRII=0;

  //FIND LIST OF UNIQUE TREES
  Unique_TreeCount_MRII=find_unique_tree_list(Unique_TreeList_MRII, Unique_FileList_MRII, &Unique_TreeCount_MRII,
  		                                        MRII_SampleName, FileNumber_MRII, snap_low, snap_high_MRII, MAX_UNIQUE_TREES);
  //Unique_TreeCount_MRII=1;
  Total_Unique_TreeCount+=Unique_TreeCount_MRII;
  //Total_Unique_TreeCount=Unique_TreeCount_MRII;
#endif


  printf("Number of Unique Trees in MR= %d\n", Unique_TreeCount_MR);
#ifdef MRII
  printf("Number of Unique Trees in MRII= %d\n", Unique_TreeCount_MRII);
#endif

	sprintf(buf, "./Samples/%s_Switch_MR_MRII_%d.dat",MR_SampleName, FileNumber_MR);
	if(!(fa = fopen(buf, "w")))
   	     {
   	       printf("can't open file %d\n", __LINE__);
   	       exit(1);
   	     }
	fprintf(fa,"%d",Unique_TreeCount_MR);
  fclose(fa);




//READ HEADERS TO constrcut unique_totNHalos and unique_TreeNHalos
  	printf("\nReading headers: load_tree_table\n\n");
    /* keep track of relevant values needed for the headers:
     * trees: tree_Ntrees(int), tree_totNHalos(int), tree_TreeNHalos(tree_Ntrees)
     * tree_dbids: nothing*/
	  Unique_TreeNHalos = malloc(sizeof(int) * Total_Unique_TreeCount);
	  Unique_TreeFirstHalo = malloc(sizeof(int) * Total_Unique_TreeCount);
	  Unique_totNHalos=0;
	  Unique_totNHalos_MR=0;
	  if(Unique_TreeCount_MR>MaxTrees)
	  	Unique_TreeCount_MR=MaxTrees;
	  for (treenr = 0; treenr < Unique_TreeCount_MR; treenr++)
	  {
	  	LastSnapShotNr=63;
	  	load_tree_table(Unique_FileList_MR[treenr], 1);
	  	Unique_totNHalos+=TreeNHalos[Unique_TreeList_MR[treenr]];
	   	Unique_totNHalos_MR+=TreeNHalos[Unique_TreeList_MR[treenr]];
	  	Unique_TreeNHalos[treenr]=TreeNHalos[Unique_TreeList_MR[treenr]];
	  	Unique_TreeFirstHalo[treenr] = TreeFirstHalo[Unique_TreeList_MR[treenr]];
	  	free(TreeFirstHalo);
	  	free(TreeNHalos);
	  }
	  printf("Finished Reading MR headers,  Unique_totNHalos=%d\n\n",Unique_totNHalos);
#ifdef MRII
	  for (treenr = 0; treenr < Unique_TreeCount_MRII; treenr++)
	  {
	  	LastSnapShotNr=67;
	  	load_tree_table(Unique_FileList_MRII[treenr], 2);
	  	Unique_totNHalos+=TreeNHalos[Unique_TreeList_MRII[treenr]];
	  	Unique_TreeNHalos[Unique_TreeCount_MR+treenr]=TreeNHalos[Unique_TreeList_MRII[treenr]];
	  	Unique_TreeFirstHalo[Unique_TreeCount_MR+treenr] = TreeFirstHalo[Unique_TreeList_MRII[treenr]];
	  	free(TreeFirstHalo);
	  	free(TreeNHalos);
	  }
		printf("\nFinished Reading headers,  Unique_totNHalos MR+MRII=%d\n\n",Unique_totNHalos);
#endif




 //READ AND WRITE TREE FILES
#ifdef TREES
		printf("\n\n doing tree files\n\n\n");

		//Open Output Files
#ifdef MR
		sprintf(buf, "%s/trees_0%d.%d%d", OutputDir, Nsnaps-1, sample_type, FileNumber_MR);
#endif
#ifdef MRII
		sprintf(buf, "%s/trees_sf1_0%d.%d%d%d", OutputDir, Nsnaps-1, sample_type,FileNumber_MR, FileNumber_MRII);
#endif
		if(!(fa = fopen(buf, "wb")))
		{
			printf("can't open file %d\n", __LINE__);
			exit(1);
		}

		fwrite(&Total_Unique_TreeCount, sizeof(int), 1, fa);
	  fwrite(&Unique_totNHalos, sizeof(int), 1, fa);
	  fwrite(Unique_TreeNHalos, sizeof(int), Total_Unique_TreeCount, fa);

	  // load and write one tree at the time (from the corresponding file)
	  for (treenr = 0; treenr < Unique_TreeCount_MR; treenr++)
	  {
	  	LastSnapShotNr=63;
	  	//to read the filenr the filenr, treeNhalos and Treefirst halo is given so we can jump into the correct
	  	//place in file Unique_FileList_MR[treenr] given by Unique_TreeFirstHalo[treenr] and read Unique_TreeNHalos[treenr] halos
	  	//Ntrees is read inside which is need to jump the header
	  	load_tree(Unique_FileList_MR[treenr], Unique_TreeNHalos[treenr], Unique_TreeFirstHalo[treenr], 1);

	  	printf("doing tree_MR for Treenr=%d (of %d), Tree Number of Halos=%d\n",treenr, Unique_TreeCount_MR, Unique_TreeNHalos[treenr]);
	  	//printf("FIRSTPROG=%d firsthalo=%d\n",Halo[0].FirstProgenitor,Halo[0].FirstHaloInFOFgroup);

	  	//write to tree file
	  	if(fwrite(Halo, sizeof(struct halo_data), Unique_TreeNHalos[treenr], fa) !=  Unique_TreeNHalos[treenr])
	  		printf("error in fwrite\n");

	  	free(Halo);

	  }
	  //and for MRII
#ifdef MRII
	  for (treenr = 0; treenr < Unique_TreeCount_MRII; treenr++)
	  {
	  	LastSnapShotNr=67;
	  	load_tree(Unique_FileList_MRII[treenr], Unique_TreeNHalos[Unique_TreeCount_MR+treenr],Unique_TreeFirstHalo[Unique_TreeCount_MR+treenr], 2 );

	  	printf("doing tree_MRII for Treenr=%d (of %d), Tree Number of Halos=%d\n",treenr, Unique_TreeCount_MRII, Unique_TreeNHalos[Unique_TreeCount_MR+treenr]);

	  	//write to tree file
	  	if(fwrite(Halo, sizeof(struct halo_data), Unique_TreeNHalos[Unique_TreeCount_MR+treenr], fa) !=  Unique_TreeNHalos[Unique_TreeCount_MR+treenr])
	  		printf("error in fwrite\n");

	  	free(Halo);
	  }
#endif

	  //close tree output file
	  fclose(fa);

	  printf("\ndone tree!\n");
#endif


#ifdef TREE_IDS
	//READ AND WRITE TREE_IDS FILES
		printf("\n\n doing tree_IDS files MR\n\n");

		//Open Output Files
#ifdef MR
		sprintf(buf,"%s/tree_dbids_0%d.%d%d", OutputDir, Nsnaps-1, sample_type, FileNumber_MR);
#endif
#ifdef MRII
		sprintf(buf, "%s/tree_sf1_dbids_0%d.%d%d%d", OutputDir, Nsnaps-1, sample_type, FileNumber_MR, FileNumber_MRII);
#endif

		if(!(fb = fopen(buf, "wb")))
		{
			printf("can't open file %d\n", __LINE__);
			exit(1);
		}

	  // load and write one tree at the time (from the corresponding file)
		for (treenr = 0; treenr < Unique_TreeCount_MR; treenr++)
	  {
			printf("Doing treeIDS_MR for Treenr=%d (of %d), Tree Number of Halos=%d\n",
					treenr, Unique_TreeCount_MR, Unique_TreeNHalos[treenr]);
			LastSnapShotNr=63;
	  	load_tree_IDs(Unique_FileList_MR[treenr], Unique_TreeNHalos[treenr], Unique_TreeFirstHalo[treenr], 1);
	  	//write to tree_dbids
	  	FOFIDs = malloc(sizeof(long long) * Unique_TreeNHalos[treenr]);
	  	int i;
	  	for(i=0;i<Unique_TreeNHalos[treenr];i++)
	  		FOFIDs[i]=HaloIDs[i].FirstHaloInFOFgroup;

	  	//currently we only write HaloIDs[i].FirstHaloInFOFgroup since we dont need anything else
	  	if(fwrite(FOFIDs, sizeof(long long), Unique_TreeNHalos[treenr], fb) !=  Unique_TreeNHalos[treenr])
	  		printf("error in fwrite\n");
	  	free(FOFIDs);
	  	//EVERYTHING
	  	//if(fwrite(HaloIDs, sizeof(struct halo_ids_data), Unique_TreeNHalos[treenr], fb) !=  Unique_TreeNHalos[treenr])
	  	//	printf("error in fwrite\n");
	  	free(HaloIDs);
	  }
//and for MRII
#ifdef MRII
		printf("\n\ doing tree_IDS files MRII\n\n");
	  for (treenr = 0; treenr < Unique_TreeCount_MRII; treenr++)
	  {
	  	printf("Doing treeIDS_MRII for Treenr=%d (of %d), Tree Number of Halos=%d\n",treenr, Unique_TreeCount_MRII, Unique_TreeNHalos[Unique_TreeCount_MR+treenr]);
	  	LastSnapShotNr=67;
	  	load_tree_IDs(Unique_FileList_MRII[treenr], Unique_TreeNHalos[Unique_TreeCount_MR+treenr], Unique_TreeFirstHalo[Unique_TreeCount_MR+treenr],2 );
	  	//write to tree_dbids - JUST FOFIDs
	  	FOFIDs = malloc(sizeof(long long) * Unique_TreeNHalos[Unique_TreeCount_MR+treenr]);
	  	int i;
	  	for(i=0;i<Unique_TreeNHalos[Unique_TreeCount_MR+treenr];i++)
	  		FOFIDs[i]=HaloIDs_MRII[i].FirstHaloInFOFgroup;

	  	if(fwrite(FOFIDs, sizeof(long long), Unique_TreeNHalos[Unique_TreeCount_MR+treenr], fb) !=  Unique_TreeNHalos[Unique_TreeCount_MR+treenr])
	  		printf("error in fwrite\n");
	  	free(FOFIDs);

	  	//EVERYTHING
	  	//if(fwrite(HaloIDs, sizeof(struct halo_ids_data), Unique_TreeNHalos[Unique_TreeCount_MR+treenr], fb) !=  Unique_TreeNHalos[Unique_TreeCount_MR+treenr])
	  	//	printf("error in fwrite\n");

	  	free(HaloIDs_MRII);
	  }
#endif

	  //close tree_IDS output file
	  fclose(fb);
	  printf("\ndone treeIds!\n\n\n");
#endif





#ifdef AUX

	  printf("Doing Aux\n");

//DO AUX FILE
	  ////Unique_totNHalos=Unique_totNHalos_MR;
	  int snap;
#ifdef MRII
	  Nsnaps=68;
	  int snap_jump=4;
#else
	  Nsnaps=64;
	  int snap_jump=0;
#endif

	  Aux_CountIDs_snap = malloc(sizeof(int) * Nsnaps);
	  Aux_OffsetIDs_snap = malloc(sizeof(int) * Nsnaps);
	  ////Aux_CountIDs_snaptree = malloc(sizeof(int) * Nsnaps * Unique_TreeCount_MR);
	  ////Aux_OffsetIDs_snaptree = malloc(sizeof(int) * Nsnaps * Unique_TreeCount_MR);
	  Aux_CountIDs_snaptree = malloc(sizeof(int) * Nsnaps * Total_Unique_TreeCount);
	  Aux_OffsetIDs_snaptree = malloc(sizeof(int) * Nsnaps * Total_Unique_TreeCount);
	  Aux_CountIDs_halo = malloc(sizeof(int) * Unique_totNHalos);
	  Aux_OffsetIDs_halo = malloc(sizeof(int) * Unique_totNHalos);

	  //printf("Nhalos=%d\n",Unique_totNHalos);

	  /*list of ids/pos/vel stored per snap, per tree, per halo
	   * offsets give the position of tree and halo in the list for that snapshot */

	  /*********************************************************************************************
	   *
	   * first just for the tree, the halos will be done later.
	   * CountIDs_halo and  OffsetIDs_halo ARE LISTED SEQUENTIALLY,
	   * NOT IN THE SNAPSHOT STRUCTURE IN WHICH THE DATA IS STORED
	   *
	   *********************************************************************************************/

	  /* CountIDs_snaptree - number of Ids for each tree at each snapshot (TotSnaps * Ntrees)
	   * OffsetIDs_snaptree - offset for each tree at each snapshot (TotSnaps * Ntrees);
	   * CountIDs_halo - Number of Ids per halo (NtotHalos)
	   * OffsetIDs_halo - global offset of each halo (int). */

	  int aux_tree_offset[Nsnaps]; // to contain, for each snapshot, the offset as we go through trees
	  int TotIds_MR;

	  for(i=0;i<Nsnaps;i++)
	  	aux_tree_offset[i]=0; // global

	  for (treenr = 0; treenr < Unique_TreeCount_MR; treenr++)
	  {
	  	//load aux data
	  	LastSnapShotNr=63;
	  	load_aux_header(Unique_FileList_MR[treenr],1);
	  	//load_all_auxdata(Unique_FileList_MR[treenr],1);

	  	int *header = TreeAuxData;
	  	TotIds_MR = header[1];
	  	AUXNtrees = header[2];
	  	TotSnaps = header[3];
	  	//if(TotSnaps != Nsnaps) printf("TotSnaps != Nsnaps !!!!\n\n");
     //printf("trees=%d snaps=%d\n",AUXNtrees,TotSnaps);

	  	printf("Initial Aux for Treenr=%d (of %d) MR NIds=%d\n",treenr, Unique_TreeCount_MR,TotIds_MR);

	  	//each array is basically made to start at a given place and have everything else bellow
	  	CountIDs_snaptree = header + 4 + 2 * TotSnaps;
	  	OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * AUXNtrees;

#ifdef MRII
	  	for(snap=0;snap<4;snap++)
	  	{
	  		Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+treenr] = 0;
	  		Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+treenr] = 0;
	  	}
#endif
	  	for(snap=snap_jump;snap<Nsnaps;snap++)
	  	{
	  		//for this snapshot get the number of ParticleIds in this tree
	  		Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+treenr]
	  		                      = CountIDs_snaptree[AUXNtrees*(snap-snap_jump) + Unique_TreeList_MR[treenr]];
	  		//for this snapshot get the offset of this tree (how many particles written before)
	  		Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+treenr] = aux_tree_offset[snap];
	  		//the offset is now being defined for each snapshot, based on the number of particles already counted
	  		aux_tree_offset[snap]+=Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap + treenr];
	  		//printf("snap=%d ofset MR=%d\n",snap, aux_tree_offset[k]);
	  	}
	  	free(TreeAuxData);
	  }

#ifdef MRII
	  int TotIds_MRII;
	  for (treenr = 0; treenr < Unique_TreeCount_MRII; treenr++)
	  {
	  	//load aux data
	  	LastSnapShotNr=67;
	  	load_aux_header(Unique_FileList_MRII[treenr],2);

	  	int *header = TreeAuxData;
	  	TotIds_MRII = header[1];
	  	AUXNtrees = header[2];
	  	TotSnaps = header[3];
	  	//if(TotSnaps != Nsnaps) printf("TotSnaps != Nsnaps !!!!\n\n");
	  	//printf("trees=%d snaps=%d\n",AUXNtrees,TotSnaps);

	  	printf("Initial Aux for Treenr=%d (of %d) MRII NIds=%d\n",treenr, Unique_TreeCount_MRII, TotIds_MRII);

	  	//each array is basically made to start at a given place and have everything else bellow
	  	CountIDs_snaptree = header + 4 + 2 * TotSnaps;
	  	OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * AUXNtrees;


	  	for(snap=0;snap<Nsnaps;snap++)
	  	{
	  		//for this snapshot get the number of ParticleIds in this tree
	  		Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr]
	  		                      = CountIDs_snaptree[AUXNtrees*snap + Unique_TreeList_MRII[treenr]];
	  		//for this snapshot get the offset of this tree (how many particles written before)
	  		Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr] = aux_tree_offset[snap];
	  		//the offset is now being defined for each snapshot, based on the number of particles already counted
	  		aux_tree_offset[snap]+=Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr];
	  	}
	  	//printf("ofset MRII=%d\n",aux_tree_offset[Nsnaps-1]);
	  	free(TreeAuxData);
	  }
#endif

	  //Now that we now how many particles all the trees have we can write how many particles are in each snapshot

	  //offset for the trees was done inside each snapshot, now do it globally
	  //for each tree the offset will be the offset due to other trees in current snap + the offset of the current snapshot
	  int snap_offset=0, snap_part_count=0;

	  for(snap=0;snap<Nsnaps;snap++)
	  {
	  	snap_part_count=0;

	  	for (treenr = 0; treenr < Total_Unique_TreeCount; treenr++)
	  	{
	  		Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+treenr]+=snap_offset;
	  		snap_part_count+= Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+treenr];
	  	}

	  	Aux_CountIDs_snap[snap] = snap_part_count;
	  	Aux_OffsetIDs_snap[snap] = snap_offset;
	  	//printf("snap=%d count=%d offset=%d\n",snap,Aux_CountIDs_snap[snap],Aux_OffsetIDs_snap[snap]);
	  	snap_offset+=snap_part_count;
	  	printf("snap_part_count=%d snap_offset=%d\n",snap_part_count,snap_offset);
	  }


	  //The total number of IDs is basically the offset after the last tree as been counted
	  Aux_TotIds=snap_offset;
	  printf("Tot IDS=%lld\n",Aux_TotIds);

	  printf("\nParticle data done for trees, now for halos and then write\n\n");

    /*************************************************************************
    ***                                                                   ***
    ***   At this moment we have the counts and offsets for the trees:    ***
    ***                                                                   ***
    ***         Aux_CountIDs_snaptree & Aux_OffsetIDs_snaptree            ***
    ***                                                                   ***
    *************************************************************************/

    //now get the halo particle counts and offsets
	  int halonr;
	  int aux_count_halos=0;
	  int *aux_count_part_snaptree;
	  long long *MyIdList;
	  float *MyPosList, *MyVelList;
	  int MyOffsetIDs;
	  //Aux_TotIds=100000000;
	  MyIdList = malloc(sizeof(long long) * Aux_TotIds);
	  MyPosList = malloc(sizeof(float) * 3 * Aux_TotIds);
	  MyVelList = malloc(sizeof(float) * 3 * Aux_TotIds);


	  for (treenr = 0; treenr < Unique_TreeCount_MR; treenr++)
	  {

	  	printf("Doing Aux tree %d (of %d)\n",treenr,Unique_TreeCount_MR);
	  	//load aux data Unique_FileList_MR[treenr]
	  	LastSnapShotNr=63;
	  	load_tree_table(Unique_FileList_MR[treenr],1);
	  	load_tree(Unique_FileList_MR[treenr], TreeNHalos[Unique_TreeList_MR[treenr]], TreeFirstHalo[Unique_TreeList_MR[treenr]],1);
	  	load_all_auxdata(Unique_FileList_MR[treenr],1);

	  	//printf("Final Aux Treenr=%d (of %d), Tree Number of Halos=%d\n",treenr, Unique_TreeCount_MR, TreeNHalos[Unique_TreeList_MR[treenr]]);

	  	int *aheader = TreeAuxData;

	  	NtotHalos = aheader[0];
	  	TotIds = aheader[1];
	  	AUXNtrees = aheader[2];
	  	TotSnaps = aheader[3];

	  	aux_count_part_snaptree = malloc(sizeof(int) * TotSnaps);

	  	//printf("NtoHalos=%d TotIds=%d Ntrees=%d Snaps=%d\n",NtotHalos,TotIds, AUXNtrees,TotSnaps);

	  	//each array is basically made to start at a given place and have everything else bellow
	  	CountIDs_snaptree = aheader + 4 + 2 * TotSnaps;
	  	OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * AUXNtrees;
	  	CountIDs_halo = OffsetIDs_snaptree + TotSnaps * AUXNtrees;
	  	OffsetIDs_halo = CountIDs_halo + NtotHalos;

	  	for(snap=0;snap<Nsnaps;snap++)
	  		aux_count_part_snaptree[snap] = 0;

	  	for(halonr=0;halonr<TreeNHalos[Unique_TreeList_MR[treenr]];halonr++)
	  	{
	  		snap=Halo[halonr].SnapNum+snap_jump;
	  		//if(Halo[halonr].MostBoundID==2743514)
	  		//	printf("tree=%d halonr=%d snap=%d found!!!\n",treenr,halonr,snap);
	  		if(CountIDs_halo[TreeFirstHalo[Unique_TreeList_MR[treenr]] + halonr] >0)
	  			Aux_OffsetIDs_halo[aux_count_halos+halonr] = aux_count_part_snaptree[snap]
	  			                                       + Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+treenr];
	  		else
	  			Aux_OffsetIDs_halo[aux_count_halos+halonr] = -1;
	  		Aux_CountIDs_halo[aux_count_halos+halonr] = CountIDs_halo[TreeFirstHalo[Unique_TreeList_MR[treenr]] + halonr];
	  		aux_count_part_snaptree[snap]+=CountIDs_halo[TreeFirstHalo[Unique_TreeList_MR[treenr]] + halonr];
	  	}

	  	for(snap=snap_jump;snap<Nsnaps;snap++)
	  	{

	  		size_t header_offset = 4 * sizeof(int) + 2 * TotSnaps * sizeof(int) + 2 * TotSnaps * AUXNtrees * sizeof(int) + 2
	  				* sizeof(int) * NtotHalos;

	  		IdList = (long long *) (TreeAuxData + header_offset);
	  		//IdList = (long long *) (TreeAuxData);
	  		PosList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long));
	  		VelList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long) + TotIds * 3 * sizeof(float));

	  		OffsetIDs = OffsetIDs_snaptree[(snap-snap_jump) * AUXNtrees + Unique_TreeList_MR[treenr]];

	  		IdList += OffsetIDs;
	  		PosList += 3 * OffsetIDs;
	  		VelList += 3 * OffsetIDs;
	  		//printf("Pos red, offsetIds_snaptree=%d Count=%d\n",Aux_OffsetIDs_snaptree[snap * Unique_TreeCount_MR +treenr],Aux_CountIDs_snaptree[snap*Unique_TreeCount_MR + treenr]);
	  		MyOffsetIDs = Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap +treenr];

	  		for(index=MyOffsetIDs; index< MyOffsetIDs+Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap + treenr];index++)
	  		{
	  			MyIdList[index]= IdList[index-MyOffsetIDs];
	  			//if(MyIdList[index]==2743514) printf("found!!!\n");
	  			for(j=0;j<3;j++)
	  			{
	  				MyPosList[3*index+j]=PosList[3*(index-MyOffsetIDs)+j];
	  				MyVelList[3*index+j]=VelList[3*(index-MyOffsetIDs)+j];
	  			}
	  		}
	  		//printf("MyPos assigned\n");
	  	}//loop on snaps
	  	//printf("index=%d\n",MyOffsetIDs);

	  	aux_count_halos+=TreeNHalos[Unique_TreeList_MR[treenr]];

	  	free(TreeFirstHalo);
	  	free(TreeNHalos);
	  	free(Halo);

	  	free(TreeAuxData);
	   	free(aux_count_part_snaptree);
	  }//loop on trees


#ifdef MRII
	  for (treenr = 0; treenr < Unique_TreeCount_MRII; treenr++)
	 	  {

	 	  	printf("Doing Aux tree %d (of %d)\n",treenr,Unique_TreeCount_MRII);
	 	  	//load aux data Unique_FileList_MR[treenr]
	 	  	LastSnapShotNr=67;
	 	  	load_tree_table(Unique_FileList_MRII[treenr],2);
	 	  	load_tree(Unique_FileList_MRII[treenr], TreeNHalos[Unique_TreeList_MRII[treenr]], TreeFirstHalo[Unique_TreeList_MRII[treenr]],2);
	 	  	load_all_auxdata(Unique_FileList_MRII[treenr],2);

	 	  	//printf("Final Aux Treenr=%d (of %d), Tree Number of Halos=%d\n",treenr, Unique_TreeCount_MR, TreeNHalos[Unique_TreeList_MR[treenr]]);

	 	  	int *aheader = TreeAuxData;

	 	  	NtotHalos = aheader[0];
	 	  	TotIds = aheader[1];
	 	  	AUXNtrees = aheader[2];
	 	  	TotSnaps = aheader[3];

	 	  	aux_count_part_snaptree = malloc(sizeof(int) * TotSnaps);

	 	  	//printf("NtoHalos=%d TotIds=%d Ntrees=%d Snaps=%d\n",NtotHalos,TotIds, AUXNtrees,TotSnaps);

	 	  	//each array is basically made to start at a given place and have everything else bellow
	 	  	CountIDs_snaptree = aheader + 4 + 2 * TotSnaps;
	 	  	OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * AUXNtrees;
	 	  	CountIDs_halo = OffsetIDs_snaptree + TotSnaps * AUXNtrees;
	 	  	OffsetIDs_halo = CountIDs_halo + NtotHalos;

	 	  	for(snap=0;snap<Nsnaps;snap++)
	 	  		aux_count_part_snaptree[snap] = 0;


	 	  	for(halonr=0;halonr<TreeNHalos[Unique_TreeList_MRII[treenr]];halonr++)
	 	  	{
	 	  		snap=Halo[halonr].SnapNum;
	 	  		//if(Halo[halonr].MostBoundID==2743514)
	 	  		//	printf("tree=%d halonr=%d snap=%d found!!!\n",treenr,halonr,snap);
	 	  		if(CountIDs_halo[TreeFirstHalo[Unique_TreeList_MRII[treenr]] + halonr] >0)
	 	  			Aux_OffsetIDs_halo[aux_count_halos+halonr] = aux_count_part_snaptree[snap] +
	 	  													Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr];
	 	  		else
	 	  			Aux_OffsetIDs_halo[aux_count_halos+halonr] = -1;
	 	  		Aux_CountIDs_halo[aux_count_halos+halonr] = CountIDs_halo[TreeFirstHalo[Unique_TreeList_MRII[treenr]] + halonr];
	 	  		aux_count_part_snaptree[snap]+=CountIDs_halo[TreeFirstHalo[Unique_TreeList_MRII[treenr]] + halonr];
	 	  	}
	 	  	for(snap=0;snap<Nsnaps;snap++)
	 	  	{

	 	  		size_t header_offset = 4 * sizeof(int) + 2 * TotSnaps * sizeof(int) + 2 * TotSnaps * AUXNtrees * sizeof(int) + 2
	 	  				* sizeof(int) * NtotHalos;

	 	  		IdList = (long long *) (TreeAuxData + header_offset);
	 	  		//IdList = (long long *) (TreeAuxData);
	 	  		PosList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long));
	 	  		VelList = (float *) (TreeAuxData + header_offset + TotIds * sizeof(long long) + TotIds * 3 * sizeof(float));

	 	  		OffsetIDs = OffsetIDs_snaptree[snap * AUXNtrees + Unique_TreeList_MRII[treenr]];

	 	  		IdList += OffsetIDs;
	 	  		PosList += 3 * OffsetIDs;
	 	  		VelList += 3 * OffsetIDs;
	 	  		//printf("Pos red, offsetIds_snaptree=%d Count=%d\n",Aux_OffsetIDs_snaptree[snap * Unique_TreeCount_MR +treenr],Aux_CountIDs_snaptree[snap*Unique_TreeCount_MR + treenr]);
	 	  		MyOffsetIDs = Aux_OffsetIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr];

	 	  		for(index=MyOffsetIDs; index< MyOffsetIDs +
	 	  												Aux_CountIDs_snaptree[Total_Unique_TreeCount*snap+Unique_TreeCount_MR+treenr];index++)
	 	  		{
	 	  			MyIdList[index]= IdList[index-MyOffsetIDs];
	 	  			//if(MyIdList[index]==2743514) printf("found!!!\n");
	 	  			for(j=0;j<3;j++)
	 	  			{
	 	  				MyPosList[3*index+j]=PosList[3*(index-MyOffsetIDs)+j];
	 	  				MyVelList[3*index+j]=VelList[3*(index-MyOffsetIDs)+j];
	 	  			}
	 	  		}

	 	  		//printf("MyPos assigned\n");
	 	  	}//loop on snaps
	 	  	//printf("index=%d\n",MyOffsetIDs);

	 	  	aux_count_halos+=TreeNHalos[Unique_TreeList_MRII[treenr]];

	 	  	free(TreeFirstHalo);
	 	  	free(TreeNHalos);
	 	  	free(Halo);

	 	  	free(TreeAuxData);

	 	  	free(aux_count_part_snaptree);
	 	  }//loop on trees

#endif

 //write aux
#ifdef MR
 sprintf(buf,"%s/treeaux_0%d.%d%d", OutputDir, Nsnaps-1,sample_type,FileNumber_MR);
#endif
#ifdef MRII
 sprintf(buf, "%s/treeaux_sf1_0%d.%d%d%d", OutputDir, Nsnaps-1,sample_type,FileNumber_MR, FileNumber_MRII);
#endif

 if(!(fc = fopen(buf, "wb")))
 {
	 printf("can't open file %d\n", __LINE__);
	 exit(1);
 }

 printf("Nhalos=%d\n",aux_count_halos);
 fwrite(&Unique_totNHalos, sizeof(int), 1, fc);
 fwrite(&Aux_TotIds, sizeof(int), 1, fc);
 fwrite(&Total_Unique_TreeCount, sizeof(int), 1, fc);
 fwrite(&Nsnaps, sizeof(int), 1, fc);

 fwrite(Aux_CountIDs_snap, sizeof(int), Nsnaps, fc);
 fwrite(Aux_OffsetIDs_snap, sizeof(int), Nsnaps, fc);
 fwrite(Aux_CountIDs_snaptree, sizeof(int), Nsnaps * Total_Unique_TreeCount, fc);
 fwrite(Aux_OffsetIDs_snaptree, sizeof(int), Nsnaps * Total_Unique_TreeCount, fc);
 fwrite(Aux_CountIDs_halo, sizeof(int), aux_count_halos, fc);
 fwrite(Aux_OffsetIDs_halo, sizeof(int), aux_count_halos, fc);

 fwrite(MyIdList, sizeof(long long), Aux_TotIds, fc);
 fwrite(MyPosList, sizeof(float), 3 * Aux_TotIds, fc);
 fwrite(MyVelList, sizeof(float), 3 * Aux_TotIds, fc);


  fclose(fc);

  printf("\nTotal Number of Particles=%d\n",Aux_TotIds);
#endif

  free(Unique_FileList_MR);
  free(Unique_TreeList_MR);
#ifdef MRII
  free(Unique_FileList_MRII);
  free(Unique_TreeList_MRII);
#endif
  free(Unique_TreeNHalos);
  free(Unique_TreeFirstHalo);


  printf("\n\ndone");

}

int find_unique_tree_list(int *Unique_TreeList, int *Unique_FileList, int *Unique_TreeCount,
		                      char *Sample_Name, int FileNumber, int snap_low, int snap_high, int MAX_UNIQUE_TREES)
{
  int i, j, treenr, snap;
  int Non_Unique_TreeCount, *Non_Unique_FileList, *Non_Unique_TreeList;
  float weight;
  long long haloid;
  FILE *fa;
  char buf[1000];

	for(snap=snap_low;snap<=snap_high;snap++)
  {
		//if(snap==53 || snap==37 || snap==30 || snap==25)
  	if(snap<70)
  	{
  	sprintf(buf, "./Samples/%s_%d%d.dat",Sample_Name, FileNumber, snap);
  	if(!(fa = fopen(buf, "r")))
  	{
  		printf("can't open file place 1 `%s'\n", buf);
  		exit(1);
  	}
    //printf("doing file %s\n",buf);
  	fscanf(fa, "%d \n", &Non_Unique_TreeCount);
  	if(Non_Unique_TreeCount>MAX_UNIQUE_TREES)
  	{
  		printf("Non_Unique_TreeCount>MAX_UNIQUE_TREES\n Increase MAX_UNIQUE_TREES\n");
  		exit(0);
  	}
  	//printf("\n\nNTrees for snap %d = %d\n", snap, Non_Unique_TreeCount);


  	Non_Unique_FileList = malloc(sizeof(int) * Non_Unique_TreeCount);
  	Non_Unique_TreeList = malloc(sizeof(int) * Non_Unique_TreeCount);


  	/* read file with Halo IDS, Tree Nr, File Nr and weights for each halo
  	 * weights are only used by the MCMC not here)*/
  	for(i=0;i<Non_Unique_TreeCount;i++)
  	{
  		fscanf(fa, "%lld %d %d %f\n", &haloid, &Non_Unique_TreeList[i], &Non_Unique_FileList[i], &weight);
  		//printf("i=%d Treelist=%d Filelist=%d\n",i,Non_Unique_TreeList[i], Non_Unique_FileList[i]);
  	}

  	fclose(fa);


 	  for (treenr = 0; treenr < Non_Unique_TreeCount; treenr++)
 	  {
 	  	if(treenr==0 && snap==snap_low)
 	  	{
 	  		*Unique_TreeCount+=1;
 	  		Unique_TreeList[0]=Non_Unique_TreeList[treenr];
 	  		Unique_FileList[0]=Non_Unique_FileList[treenr];

 	  	}
 	  	else
 	  		for(j=0;j<*Unique_TreeCount;j++)
 	  		{
 	  			//printf("treenr=%d j=%d tree=%d Unique=%d\n",treenr,j, TreeList[treenr],Unique_TreeList[j]);

 	  			// If the tree has already been done, jump to next tree in the list
 	  			if(treenr==Non_Unique_TreeCount ||
 	  					(Non_Unique_TreeList[treenr]==Unique_TreeList[j] && Non_Unique_FileList[treenr]==Unique_FileList[j]))
 	  				break;
 	  			else
 	  				if(j==*Unique_TreeCount-1 && Non_Unique_TreeList[treenr]!=Unique_TreeList[*Unique_TreeCount-1]
 	  				                          && Non_Unique_FileList[treenr]!=Unique_FileList[*Unique_TreeCount-1])
 	  				{
 	  					*Unique_TreeCount+=1;
 	  					Unique_TreeList[*Unique_TreeCount-1]=Non_Unique_TreeList[treenr];
 	  					Unique_FileList[*Unique_TreeCount-1]=Non_Unique_FileList[treenr];
 	  					//printf("New Unique tree=%d in position %d\n",Unique_TreeList[Unique_TreeCount-1],Unique_TreeCount-1);
 	  					break;
 	  				}
 	  		}

 	  	if(treenr==Non_Unique_TreeCount)
 	  		break;

 	  }//end loop on non unique trees

  	}//end if snap==
  }//end loop on snap files
  return *Unique_TreeCount;
}


