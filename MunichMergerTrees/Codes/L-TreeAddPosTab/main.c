#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

int SnapNum;

struct bigtab_data
{
  int index;
  long long ID;
  float Pos[3];
  float Vel[3];
} *BigTab;

int swap_file=8;

char label[4];
int  nextblock;

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data,int n,int m)
{
  int i,j;
  char old_data[16];

  if(swap_file!=8)
    {
      for(j=0;j<n;j++)
	{
          memcpy(&old_data[0],&data[j*m],m);
          for(i=0;i<m;i++)
            {
              data[j*m+i]=old_data[m-i-1];
	    }
	}
    }
}

/*------------------------------------------------------------------*/
/*----------- procedure to swap header if needed -------------------*/
/*------------------------------------------------------------------*/

void swap_header(void)
{
   swap_Nbyte((char*)&header.npart,6,4);
   swap_Nbyte((char*)&header.mass,6,8);
   swap_Nbyte((char*)&header.time,1,8);
   swap_Nbyte((char*)&header.redshift,1,8);
   swap_Nbyte((char*)&header.flag_sfr,1,4);
   swap_Nbyte((char*)&header.flag_feedback,1,4); 
   swap_Nbyte((char*)&header.npartTotal,6,4);
   swap_Nbyte((char*)&header.flag_cooling,1,4);
   swap_Nbyte((char*)&header.num_files,1,4);  
   swap_Nbyte((char*)&header.BoxSize,1,8);
   swap_Nbyte((char*)&header.Omega0,1,8);
   swap_Nbyte((char*)&header.OmegaLambda,1,8);
   swap_Nbyte((char*)&header.HubbleParam,1,8);
   swap_Nbyte((char*)&header.flag_stellarage,1,4);
   swap_Nbyte((char*)&header.flag_metals,1,4);
   /*
   swap_Nbyte((char*)&header.npartTotalHighWord,6,4);
   swap_Nbyte((char*)&header.flag_entropy_instead_u,1,4);
   */
}

/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void find_block(char *label,FILE *fd)
{
  int blocksize=0,blksize;
  char blocklabel[5]={"    "};

#define FBSKIP  {fread(&blksize,sizeof(int),1,fd);}

  rewind(fd);

  while(!feof(fd) && blocksize == 0)
    {
      FBSKIP;
      swap_file=blksize;
      swap_Nbyte((char*)&blksize,1,4);
      if(blksize != 8)
	{
	  printf("Incorrect Format (blksize=%d)!\n",blksize);
	  exit(1891);
	}
       else
         {
           fread(blocklabel, 4*sizeof(char), 1, fd);
           fread(&blocksize, sizeof(int), 1, fd);
           swap_Nbyte((char*)&blocksize,1,4);
	   printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
                   label[0],label[1],label[2],label[3],blocklabel,blocksize);
           FBSKIP;
           if(strncmp(label,blocklabel,4)!=0)
             { 
	       fseek(fd,blocksize,1);
	       blocksize=0;
             }
         }
    }
  if(feof(fd))
    {
       printf("Block '%c%c%c%c' not found !\n",
	      label[0],label[1],label[2],label[3]);
    }
}

int test_block(char *label,FILE *fd)
{
  int blocksize=0,blksize;
  char blocklabel[5]={"    "};
  size_t len;

#define TBSKIP  {len = fread(&blksize,sizeof(int),1,fd);}

  rewind(fd);

  while(!feof(fd) && blocksize == 0)
    {
      TBSKIP;
      if(len != 1)
	{
	  printf("Skip failed ...\n");
          break;
	} 
      swap_file=blksize;
      swap_Nbyte((char*)&blksize,1,4);
      if(blksize != 8)
        {
          printf("Incorrect Format (blksize=%d)!\n",blksize);
          exit(1891);
        }
       else
         {
           fread(blocklabel, 4*sizeof(char), 1, fd);
           fread(&blocksize, sizeof(int), 1, fd);
           swap_Nbyte((char*)&blocksize,1,4);
           printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
                   label[0],label[1],label[2],label[3],blocklabel,blocksize);
           TBSKIP;
           if(strncmp(label,blocklabel,4)!=0)
             {
               printf("seek: %d\n",fseek(fd,blocksize,1));
                blocksize=0;
             }
         }
    }
  if(feof(fd))
      return(-1);
  else
      return(0);
}


int main(int argc, char **argv)
{
  int snapfilenr, matches;

  if(argc != 3)
    {
      printf("\n  usage: L-TreeAddIDTab <parameterfile>  <snapnum> \n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);

  SnapNum = atoi(argv[2]);



  count_ids(SnapNum);

  BigTab = mymalloc(sizeof(struct bigtab_data) * TotIDs);

  load_ids(SnapNum);

  qsort(BigTab, TotIDs, sizeof(struct bigtab_data), compare_bigtab_id);


  for(snapfilenr = 0; snapfilenr < FilesPerSnapshot; snapfilenr++)
    {
      load_pos_data(SnapNum, snapfilenr);

      printf("%s_%03d.%d  .....  ", SnapshotFileBase, SnapNum, snapfilenr);
      fflush(stdout);

      matches = do_match(SnapNum);

      free_pos_data();
    }

  qsort(BigTab, TotIDs, sizeof(struct bigtab_data), compare_bigtab_index);

  save_ids(SnapNum);

  myfree(BigTab);

  return 0;
}


void count_ids(int snapnum)
{
  int filenr, totHalos, totIds, totTrees, totSnaps;
  char buf[1000];
  int *countIDs_snap;
  FILE *fd;

  for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
      sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("[c] can't open file `%s'\n", buf);
	  exit(1);
	}

      fread(&totHalos, 1, sizeof(int), fd);
      fread(&totIds, 1, sizeof(int), fd);
      fread(&totTrees, 1, sizeof(int), fd);
      fread(&totSnaps, 1, sizeof(int), fd);

      countIDs_snap = mymalloc(sizeof(int) * totSnaps);
      fread(countIDs_snap, totSnaps, sizeof(int), fd);

      TotIDs += countIDs_snap[snapnum];

      myfree(countIDs_snap);

      fclose(fd);
    }

  printf("TotIDs = %d\n", TotIDs);
}




void load_ids(int snapnum)
{
  int i, filenr, totHalos, totIds, totTrees, totSnaps;
  int count;
  char buf[1000];
  int *countIDs_snap, *offsetIDs_snap;
  long long *ids;
  FILE *fd;


  count = 0;

  for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
      sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("[c] can't open file `%s'\n", buf);
	  exit(1);
	}
      
      printf("load %s\n", buf);

      fread(&totHalos, 1, sizeof(int), fd);
      fread(&totIds, 1, sizeof(int), fd);
      fread(&totTrees, 1, sizeof(int), fd);
      fread(&totSnaps, 1, sizeof(int), fd);

      countIDs_snap = mymalloc(sizeof(int) * totSnaps);
      offsetIDs_snap = mymalloc(sizeof(int) * totSnaps);
      fread(countIDs_snap, totSnaps, sizeof(int), fd);
      fread(offsetIDs_snap, totSnaps, sizeof(int), fd);


      fseek(fd, 2 * sizeof(int) * totSnaps * totTrees + 2 * sizeof(int) * totHalos +
	    offsetIDs_snap[snapnum] * sizeof(long long), SEEK_CUR);

      ids = mymalloc(countIDs_snap[snapnum] * sizeof(long long));

      fread(ids, sizeof(long long), countIDs_snap[snapnum], fd);

      for(i = 0; i < countIDs_snap[snapnum]; i++)
	{
	  BigTab[count].ID = ids[i];
	  BigTab[count].index = count;
	  count++;
	}

      myfree(ids);
      myfree(offsetIDs_snap);
      myfree(countIDs_snap);

      fclose(fd);
    }
}


void save_ids(int snapnum)
{
  int i, k, count, countpos, filenr, totHalos, totIds, totTrees, totSnaps;
  char buf[1000];
  int *countIDs_snap, *offsetIDs_snap;
  float *pos;
  FILE *fd;


  count = 0;

  printf("saving...\n");

  for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
      sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
      if(!(fd = fopen(buf, "r+")))
	{
	  printf("[c] can't open file `%s'\n", buf);
	  exit(1);
	}

      fread(&totHalos, 1, sizeof(int), fd);
      fread(&totIds, 1, sizeof(int), fd);
      fread(&totTrees, 1, sizeof(int), fd);
      fread(&totSnaps, 1, sizeof(int), fd);

      countIDs_snap = mymalloc(sizeof(int) * totSnaps);
      offsetIDs_snap = mymalloc(sizeof(int) * totSnaps);
      fread(countIDs_snap, totSnaps, sizeof(int), fd);
      fread(offsetIDs_snap, totSnaps, sizeof(int), fd);


      fseek(fd, 2 * sizeof(int) * totSnaps * totTrees + 2 * sizeof(int) * totHalos +
	    totIds * sizeof(long long) + offsetIDs_snap[snapnum] * 3 * sizeof(float), SEEK_CUR);

      /* write coordinates */

      pos = mymalloc(countIDs_snap[snapnum] * 3 * sizeof(float));

      for(i = 0, countpos= count; i < countIDs_snap[snapnum]; i++)
	{
	  for(k = 0; k < 3; k++)
	    pos[i * 3 + k] = BigTab[countpos].Pos[k];
	  countpos++;
	}

      fwrite(pos, 3 * sizeof(float), countIDs_snap[snapnum], fd);

      /* write velocities */

      for(i = 0; i < countIDs_snap[snapnum]; i++)
	{
	  for(k = 0; k < 3; k++)
	    pos[i * 3 + k] = BigTab[count].Vel[k];
	  count++;
	}

      fseek(fd, (totIds - countIDs_snap[snapnum])* 3 * sizeof(float), SEEK_CUR);
      fwrite(pos, 3 * sizeof(float), countIDs_snap[snapnum], fd);


      myfree(pos);
      myfree(offsetIDs_snap);
      myfree(countIDs_snap);

      fclose(fd);
    }

  printf("done.\n");

}







void load_pos_data(int num, int filenr)
{
  int i, j, dummy;
  char buf[1000];
  FILE *fd;
  float *pos;
#ifdef LONGIDS
  long long *id;
#else
  int *id;
#endif
#ifdef HAVE_HDF5
  hid_t file_id;
#endif

  if(SnapFormat != 3)
    {
      if(FilesPerSnapshot > 1)
	//sprintf(buf, "%s/%s_%03d.%d", SimulationDir, SnapshotFileBase, num, filenr);
	sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", SimulationDir, num, SnapshotFileBase, num, filenr);
      else
	//sprintf(buf, "%s/%s_%03d", SimulationDir, SnapshotFileBase, num);
	sprintf(buf, "%s/snapdir_%03d/%s_%03d", SimulationDir, num, SnapshotFileBase, num);
    }
  else
    {
      if(FilesPerSnapshot > 1)
	//sprintf(buf, "%s/%s_%03d.%d", SimulationDir, SnapshotFileBase, num, filenr);
	sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5", SimulationDir, num, SnapshotFileBase, num, filenr);
      else
	//sprintf(buf, "%s/%s_%03d", SimulationDir, SnapshotFileBase, num);
	sprintf(buf, "%s/snapdir_%03d/%s_%03d.hdf5", SimulationDir, num, SnapshotFileBase, num);
    }

  if(SnapFormat != 3)
    {
      if(!(fd = fopen(buf, "r")))
	{
	  printf("[a] can't open file `%s'\n", buf);
	  exit(1);
	}
    }
  else
    {
#ifdef HAVE_HDF5
      if((file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
	{
	  printf("[a] can't open file `%s'\n", buf);
	  exit(1);
	}
#else
      printf("SnapFormat=3 but code was compiled without HDF5\n");
      exit(1);
#endif
    }

  if(SnapFormat != 3)
    {
      if(SnapFormat == 2)
	{
	  find_block("HEAD",fd);
	}
      fread(&dummy, sizeof(int), 1, fd);

      if(SnapFormat == 1)
	{
	  if (dummy!=256) swap_file = dummy; else swap_file = 8;
	}
      fread(&header, sizeof(header), 1, fd);
      swap_header();
      fread(&dummy, sizeof(int), 1, fd);
    }
  else
    {
#ifdef HAVE_HDF5
      read_hdf5_header(file_id, &header);
#else
      printf("SnapFormat=3 but code was compiled without HDF5\n");
      exit(1);
#endif
    }

  NumPart = header.npart[0] + header.npart[1] + header.npart[2] + 
            header.npart[3] + header.npart[4] + header.npart[5];

  printf("NumPart = %d\n",NumPart);

  P = mymalloc(NumPart * sizeof(struct particle_data));

  pos = mymalloc(NumPart * 3 * sizeof(float));

  if(SnapFormat != 3)
    {
      if(SnapFormat == 2)
	{
	  find_block("POS ",fd);
	}
      fread(&dummy, sizeof(int), 1, fd);

      fread(pos, 3 * sizeof(float), NumPart, fd);
      swap_Nbyte((char*)pos,3 * NumPart,4);
      fread(&dummy, sizeof(int), 1, fd);
    }
  else
    {
#ifdef HAVE_HDF5
      /* Read positions */
      read_all_types(file_id, &header, "Coordinates", H5T_NATIVE_FLOAT, pos);      
#else
      printf("SnapFormat=3 but code was compiled without HDF5\n");
      exit(1);
#endif
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Pos[j] = pos[i * 3 + j];

  if(SnapFormat != 3)
    {
      if(SnapFormat == 2)
	{
	  find_block("VEL ",fd);
	}
      fread(&dummy, sizeof(int), 1, fd);

      fread(pos, 3 * sizeof(float), NumPart, fd);  /* load velocities */
      swap_Nbyte((char*)pos,3 * NumPart,4);
      fread(&dummy, sizeof(int), 1, fd);
    }
  else
    {
#ifdef HAVE_HDF5
      /* Read velocities */
      read_all_types(file_id, &header, "Velocities", H5T_NATIVE_FLOAT, pos);
#else
      printf("SnapFormat=3 but code was compiled without HDF5\n");
      exit(1);
#endif
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Vel[j] = pos[i * 3 + j];

  myfree(pos);




#ifdef LONGIDS
  id = mymalloc(NumPart * sizeof(long long));
#else
  id = mymalloc(NumPart * sizeof(int));
#endif

  if(SnapFormat != 3)
    {
      if(SnapFormat == 2)
	{
#ifdef PATCH_IDS
	  if(test_block("IDU ",fd) == -1) find_block("ID  ",fd);
#else
	  find_block("ID  ",fd);
#endif
	}
      fread(&dummy, sizeof(int), 1, fd);
#ifdef LONGIDS
      fread(id, sizeof(long long), NumPart, fd);
      swap_Nbyte((char*)id,NumPart,8);
#else
      fread(id, sizeof(int), NumPart, fd);
      swap_Nbyte((char*)id,NumPart,4);
#endif
      fread(&dummy, sizeof(int), 1, fd);
    }
  else
    {
#ifdef HAVE_HDF5
      /* Read IDs */
      read_all_types(file_id, &header, "ParticleIDs", HDF5_ID_TYPE, id);
#else
      printf("SnapFormat=3 but code was compiled without HDF5\n");
      exit(1);
#endif
    }


  for(i = 0; i < NumPart; i++)
    P[i].ID = id[i];

  myfree(id);

  if(SnapFormat != 3)
    fclose(fd);
#ifdef HAVE_HDF5
  else
    H5Fclose(file_id);
#endif

  qsort(P, NumPart, sizeof(struct particle_data), compare_snap_id);
}



void free_pos_data(void)
{
  myfree(P);
}







int do_match(int snapnum)
{
  int i, j, k;
  int matches = 0;

  for(i = 0, j = 0; i < TotIDs; i++)
    {
      while(P[j].ID < BigTab[i].ID)
	{
	  j++;

	  if(j >= NumPart)
	    break;
	}

      if(j >= NumPart)
	break;

      if(P[j].ID == BigTab[i].ID)
	{
	  for(k = 0; k < 3; k++)
            {
              BigTab[i].Pos[k] = P[j].Pos[k];
              BigTab[i].Vel[k] = P[j].Vel[k];
            }

	  matches++;
	}
    }

  printf("matches = %d\n", matches);

  return matches;
}






int compare_snap_id(const void *a, const void *b)
{
  if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
    return -1;

  if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
    return +1;

  return 0;
}



int compare_bigtab_id(const void *a, const void *b)
{
  if(((struct bigtab_data *) a)->ID < (((struct bigtab_data *) b)->ID))
    return -1;

  if(((struct bigtab_data *) a)->ID > (((struct bigtab_data *) b)->ID))
    return +1;

  return 0;
}



int compare_bigtab_index(const void *a, const void *b)
{
  if(((struct bigtab_data *) a)->index < (((struct bigtab_data *) b)->index))
    return -1;

  if(((struct bigtab_data *) a)->index > (((struct bigtab_data *) b)->index))
    return +1;

  return 0;
}
