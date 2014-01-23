#include "output_mt.h"

void create_subfind_substruct(m_halo_wrapper_t* haloB);
void internalaux_read(clgal_aux_data_wrapper_t *aux_data, char* outputfolder);


int compare_clgal_aux_data_t_by_globalRefID(const void *v1, const void *v2);
uint64_t search_clgal_aux_data_t_for_globalRefID( uint64_t searchID, uint64_t n_array ,const void *Array );

int refdomain;
typedef struct full_tree
{
  hid_t intreeid,globalRefID;
} full_tree_t;

typedef struct id_component
{
  int snapid;
  int domainid;
  hid_t localid;
} id_component_t;

id_component_t extract_id_component(hid_t hid);

void treecrawler(hid_t hid, clgal_aux_data_wrapper_t **aux_data, int treenr, full_tree_t **fulltree, hid_t *nHalosinTree);
void complete_clgal_aux(hid_t hid, hid_t refid, clgal_aux_data_wrapper_t **aux_data, char* outputfolder);
clgal_aux_data_t* clgal_aux_data_pointer_from_globalRefID(hid_t hid, clgal_aux_data_wrapper_t **aux_data);


/* For God's sake, I need this function to point the pointer to an element of aux_data when I specify a globalRefID - Boyd */
clgal_aux_data_t* clgal_aux_data_pointer_from_globalRefID(hid_t hid, clgal_aux_data_wrapper_t **aux_data)
{
  id_component_t local_snap_data;
  clgal_aux_data_t* data;
  local_snap_data = extract_id_component(hid);
  data = &(aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid]);
  return data;
}


void generate_lgal_output(char* outputfolder, int localdomain,float *snaplist, int nSnaps, int totaldomains)
{
  clgal_aux_data_wrapper_t **aux_data;
  int i,j,itree,total_trees;
  hid_t ihalo,curid,cur_fof_id;
  hid_t *nHalosinTree;
  full_tree_t **fulltree;
  id_component_t local_snap_data;
  clgal_aux_data_t* cur_aux_data;

  /* Set up snapshot info */
  aux_data = malloc(nSnaps*sizeof(clgal_aux_data_wrapper_t*));
  refdomain = localdomain;
  for(i=1;i<nSnaps;i++)
    {
      aux_data[i] = malloc(totaldomains*sizeof(clgal_aux_data_wrapper_t));
      for(j=0;j<totaldomains;j++)
	{
	  aux_data[i][j].already_read = 0;
	  aux_data[i][j].redshift = snaplist[i];
	  aux_data[i][j].snapid = i;
	  aux_data[i][j].domainid = j;	  
	}
    }

  /* Read in last snapshot data */
  internalaux_read(&(aux_data[nSnaps-1][localdomain]), outputfolder);
  
  /* Prepare other snapshots */
  for(ihalo=0;ihalo<aux_data[nSnaps-1][localdomain].nHalos;ihalo++)
    {
      complete_clgal_aux(aux_data[nSnaps-1][localdomain].lgal_aux_halos[ihalo].globalRefID, NULLPOINT, aux_data, outputfolder);
    }

  
  /* Set up merger_tree - from root level */
  fulltree = malloc(aux_data[nSnaps-1][localdomain].nHalos*sizeof(full_tree_t *));
  nHalosinTree = calloc(aux_data[nSnaps-1][localdomain].nHalos,sizeof(hid_t));

  for(ihalo=0;ihalo<aux_data[nSnaps-1][localdomain].nHalos;ihalo++)
    {
      fulltree[ihalo] = malloc(0);
      treecrawler(aux_data[nSnaps-1][localdomain].lgal_aux_halos[ihalo].globalRefID, aux_data, (int)ihalo, fulltree, nHalosinTree);
    }


  /* Group trees into bushes */
  total_trees = aux_data[nSnaps-1][localdomain].nHalos;
  for(itree=0;itree<total_trees;itree++)
    {
      for(ihalo=0;ihalo<nHalosinTree[itree];ihalo++)
	{
	  curid = fulltree[itree][ihalo].globalRefID;
	  cur_aux_data = clgal_aux_data_pointer_from_globalRefID(curid,aux_data);
	  cur_fof_id = cur_aux_data->FirstFOF;
	  while(cur_fof_id < NULLPOINT)
	    {
	      //printf("checking halo %llu=>%llu\n",curid,cur_fof_id);
	      cur_aux_data = clgal_aux_data_pointer_from_globalRefID(cur_fof_id,aux_data);
	      if(cur_aux_data->TreeNr != itree)
		{
		  printf("checking halo %llu=>%llu\n",curid,cur_fof_id);
		  printf("moving %d => %d\n",cur_aux_data->TreeNr,itree);
		}
	      curid = cur_fof_id;
	      cur_fof_id = cur_aux_data->NextFOF;
	    }
	  //printf("checking halo %llu=>%llu\n",curid,cur_fof_id);
	}
    }


  for(ihalo=0;ihalo<aux_data[nSnaps-1][localdomain].nHalos;ihalo++)
    {
      free(fulltree[ihalo]);
    }
  free(fulltree);
  free(nHalosinTree);

  for(i=1;i<nSnaps;i++)
    {
      for(j=0;j<totaldomains;j++)
	{
	  if(aux_data[i][j].already_read == 1)
	    {
	      for(ihalo=0;ihalo<aux_data[i][j].nHalos;ihalo++)
		{
		  free(aux_data[i][j].lgal_aux_halos[ihalo].proglist);
		}
	    }  
	}
      free(aux_data[i]);     
    }
  free(aux_data);
}

/* fill in FirstProgenitors,NextProgenitor,Descendant */
void complete_clgal_aux(hid_t hid, hid_t refid, clgal_aux_data_wrapper_t **aux_data, char* outputfolder)
{
  hid_t ihalo,whalo;
  int snapid,domainid;
  hid_t localid;
  hid_t progid,nextprog,previd,curid;
  int i;
  id_component_t local_snap_data;
  snapid = hid/(uint64_t)pow(10,15);
  domainid = (hid%(uint64_t)pow(10,15))/(uint64_t)pow(10,10);
  localid = hid%(uint64_t)pow(10,10)-1;


  /* make sure to read the current catalogue */
  internalaux_read(&(aux_data[snapid][domainid]),outputfolder);

  if(aux_data[snapid][domainid].lgal_aux_halos[localid].doneaux == 1)
    {
      printf("%llu duplicated in complete aux\n",aux_data[snapid][domainid].lgal_aux_halos[localid].globalRefID);
      printf("Desc: %llu\n",aux_data[snapid][domainid].lgal_aux_halos[localid].prevDesc);
      printf("RefID : %llu\n",refid);
      
    }


  for(i=0;i<aux_data[snapid][domainid].lgal_aux_halos[localid].nprogs;i++)
    {
      if(i==0)
	{
	  curid = aux_data[snapid][domainid].lgal_aux_halos[localid].proglist[i];
	  complete_clgal_aux(curid, hid, aux_data, outputfolder);
	  local_snap_data = extract_id_component(curid);
	  internalaux_read(&(aux_data[local_snap_data.snapid][local_snap_data.domainid]), outputfolder);
	  aux_data[snapid][domainid].lgal_aux_halos[localid].FirstProgenitor = curid;
	  aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid].Descendant = hid;
	  aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid].prevDesc = hid;
	  previd = curid;
	}
      else
	{
	  curid = aux_data[snapid][domainid].lgal_aux_halos[localid].proglist[i];
	  complete_clgal_aux(curid, hid,aux_data, outputfolder);
	  local_snap_data = extract_id_component(curid);
	  internalaux_read(&(aux_data[local_snap_data.snapid][local_snap_data.domainid]), outputfolder);
	  local_snap_data = extract_id_component(previd);
	  aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid].NextProgenitor = curid;
	  local_snap_data = extract_id_component(curid);
	  aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid].Descendant = hid;
	  aux_data[local_snap_data.snapid][local_snap_data.domainid].lgal_aux_halos[local_snap_data.localid].prevDesc = hid;
	  previd = curid;
	}
    }
  aux_data[snapid][domainid].lgal_aux_halos[localid].doneaux = 1;
}

id_component_t extract_id_component(hid_t hid)
{
  id_component_t extract;
  extract.snapid = hid/(uint64_t)pow(10,15);
  extract.domainid = (hid%(uint64_t)pow(10,15))/(uint64_t)pow(10,10);
  extract.localid = hid%(uint64_t)pow(10,10)-1;
  return extract;
}

void treecrawler(hid_t hid, clgal_aux_data_wrapper_t **aux_data, int treenr, full_tree_t **fulltree, hid_t *nHalosinTree)
{
  int snapid,domainid;
  hid_t localid;
  hid_t progid,nextprog;
  snapid = hid/(uint64_t)pow(10,15);
  domainid = (hid%(uint64_t)pow(10,15))/(uint64_t)pow(10,10);
  localid = hid%(uint64_t)pow(10,10)-1;

  if(aux_data[snapid][domainid].lgal_aux_halos[localid].donetree == 1)
    {
      printf("%llu duplicated in treecrawler\n",aux_data[snapid][domainid].lgal_aux_halos[localid].globalRefID);
    }
  //printf("hid = %llu\n",hid);
  //internalaux_read(&(aux_data[snapid][domainid]), outputfolder,outputfolder);
  /* if(aux_data[snapid][domainid].lgal_aux_halos[localid].TreeNr > -1) */
  /*   { */
  /*     printf("id: %llu duplicated : original %d\n",aux_data[snapid][domainid].lgal_aux_halos[localid].globalRefID,aux_data[snapid][domainid].lgal_aux_halos[localid].TreeNr); */
  /*   } */
  if(aux_data[snapid][domainid].lgal_aux_halos[localid].TreeNr > -1)
    {
      printf("%llu %d -> %d | < %llu ... %llu\n",aux_data[snapid][domainid].lgal_aux_halos[localid].globalRefID,aux_data[snapid][domainid].lgal_aux_halos[localid].TreeNr,treenr,aux_data[snapid][domainid].lgal_aux_halos[localid].Descendant,aux_data[38][refdomain].lgal_aux_halos[treenr].globalRefID);
    }
  aux_data[snapid][domainid].lgal_aux_halos[localid].TreeNr = treenr;
  aux_data[snapid][domainid].lgal_aux_halos[localid].hidTree = nHalosinTree[treenr];


  /* add hid to fulltree */
  nHalosinTree[treenr]++;
  fulltree[treenr] = realloc(fulltree[treenr],nHalosinTree[treenr]*sizeof(full_tree_t));
  fulltree[treenr][nHalosinTree[treenr]-1].globalRefID = hid;

  progid = aux_data[snapid][domainid].lgal_aux_halos[localid].FirstProgenitor;
  //printf("firstprog: %llu => %llu\n",hid,progid);
  if(progid < NULLPOINT)
    {
      treecrawler(progid, aux_data, treenr, fulltree, nHalosinTree);
    }
  nextprog = aux_data[snapid][domainid].lgal_aux_halos[localid].NextProgenitor;
  //printf("nextprog: %llu => %llu\n",hid,nextprog);
  if(nextprog < NULLPOINT)
    {
      treecrawler(nextprog, aux_data, treenr, fulltree, nHalosinTree);
    }
  aux_data[snapid][domainid].lgal_aux_halos[localid].donetree = 1;
}

void sussexbigrun_dm_outputs( m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  sprintf(foldername,"%s/%3.3f",outputfolder,haloB->redshift);
  sprintf(command,"mkdir -p %s", foldername);
  system(command);
  sprintf(filename,"%s/%3.3f/dmdt_%d.dat",outputfolder,haloB->redshift,domainid);
  sprintf(command,"rm -f %s",filename);
  system(command);
 
  fp = fopen(filename, "w+");
  for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
    {
      if(haloB->mhalos[ihalo].used == 1)
	{
	  fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%llu\t%d\n",
		  haloB->mhalos[ihalo].Xc,
		  haloB->mhalos[ihalo].Yc,
		  haloB->mhalos[ihalo].Zc,
		  haloB->mhalos[ihalo].Mvir,
		  haloB->mhalos[ihalo].Rvir,
		  haloB->mhalos[ihalo].dm_dt,
		  haloB->mhalos[ihalo].oriID,
		  haloB->mhalos[ihalo].domainID);
	}
    }
  fclose(fp);
}

/* Need haloB->mhalos[x].ID = x */
void create_subfind_substruct(m_halo_wrapper_t* haloB)
{
  hid_t ihalo;
  hid_t upperhost,hosthalo,ref,cursub,curid;

  maphalo_to_host_mt(haloB);
  //printf("Generate Subfind structure\n");

  for(ihalo=0;ihalo<haloB->nHalos;ihalo++)
    {
      hosthalo = haloB->mhalos[ihalo].host_halo;
      upperhost = haloB->mhalos[ihalo].host_halo;
      ref = upperhost;
      while(upperhost < NULLPOINT)
	{
	  ref = upperhost;
	  upperhost = haloB->mhalos[upperhost].host_halo;
	  //printf("ref = %llu\n",ref);
	}
      if(hosthalo != ref)
	{
	  hosthalo = ref;
	}
      haloB->mhalos[ihalo].UpHalo = hosthalo;
      //printf("1 %llu hosthalo = %llu\n",ihalo,hosthalo);
      if(haloB->mhalos[ihalo].UpHalo < NULLPOINT)
	{
	  cursub = haloB->mhalos[ihalo].UpHalo;
	  while(cursub < NULLPOINT)
	    {
	      curid = haloB->mhalos[cursub].ID;
	      cursub = haloB->mhalos[cursub].NextHalo;
	    }
	  haloB->mhalos[curid].NextHalo = ihalo;
	}
      //printf("2 %llu hosthalo = %llu\n",ihalo,hosthalo);
    }
  //printf("finish linking struct\n");
  /* rename everything to globalRefID */
  for(ihalo=0;ihalo< haloB->nHalos; ihalo++)
    {
      if(haloB->mhalos[ihalo].NextHalo < NULLPOINT)
	haloB->mhalos[ihalo].NextHalo = haloB->mhalos[haloB->mhalos[ihalo].NextHalo].globalRefID;
      if(haloB->mhalos[ihalo].UpHalo < NULLPOINT)
	haloB->mhalos[ihalo].UpHalo = haloB->mhalos[haloB->mhalos[ihalo].UpHalo].globalRefID;	
    }
  //printf("finish relabel struct\n");
} 

void internalaux_outputs(m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo,whalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  int ih,len;
  float M200,Pos[3],Vel[3],VelDisp,Vmax,Spin[3];
  long long MostBoundID;
  hid_t nHalos;

 
  create_subfind_substruct(haloB);

  sprintf(foldername,"%s/%3.3f",outputfolder,haloB->redshift);
  sprintf(command,"mkdir -p %s", foldername);
  system(command);
  sprintf(filename,"%s/%3.3f/mtaux_%d.dat_bin",outputfolder,haloB->redshift,domainid);
  sprintf(command,"rm -f %s",filename);
  system(command);

  fp = fopen(filename, "wb+");
  if(fp != NULL)
    {
      /* write total halo number */
      nHalos = 0;
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      nHalos++;
	    }
	}
      fwrite(&(nHalos),sizeof(hid_t),1,fp);
      /* write nprogs */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{	  
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      fwrite(&(haloB->mhalos[ihalo].nprogs),sizeof(uint32_t),1,fp);
	    }
	}
      /* write proglist */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  //printf("Just write in : %llu --- %d\n",haloB->mhalos[ihalo].globalRefID,haloB->mhalos[ihalo].nprogs);
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      for(whalo=0; whalo < haloB->mhalos[ihalo].nprogs; whalo++)
		{	      
		  /* for(whalo=0;whalo<haloB->mhalos[ihalo].nprogs;whalo++) */
		  /* 	printf("%llu: prog  ----> %llu\n",whalo,haloB->mhalos[ihalo].proglist[whalo]); */
		  
		  fwrite(&(haloB->mhalos[ihalo].proglist[whalo]),sizeof(hid_t),1,fp);
		}
	    }
	}


      /* write globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      //printf("globalRefID: %llu\n",haloB->mhalos[ihalo].globalRefID);
	      fwrite(&(haloB->mhalos[ihalo].globalRefID),sizeof(hid_t),1,fp);
	    }
	}
      /* write FirstFOF globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      //printf("FirstFOF: %llu\n",haloB->mhalos[ihalo].UpHalo);
	      fwrite(&(haloB->mhalos[ihalo].UpHalo),sizeof(hid_t),1,fp);
	    }
	}  
      /* write NextFOF globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      //printf("NextFOF: %llu\n",haloB->mhalos[ihalo].NextHalo);
	      fwrite(&(haloB->mhalos[ihalo].NextHalo),sizeof(hid_t),1,fp);
	    }
	}  
      /* write M200 */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      M200 = haloB->mhalos[ihalo].Mvir * Msun2Gadget;
	      //printf("M200: %f\n",M200);
	      fwrite(&(M200),sizeof(float),1,fp);
	    }
	}
      /* write len */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      len = haloB->mhalos[ihalo].npart;
	      fwrite(&(len),sizeof(int),1,fp);
	    }
	}
      /* write Pos[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      Pos[0] = haloB->mhalos[ihalo].Xc*kpc2Gadget;
	      Pos[1] = haloB->mhalos[ihalo].Yc*kpc2Gadget;
	      Pos[2] = haloB->mhalos[ihalo].Zc*kpc2Gadget;
	      fwrite(Pos,sizeof(float),3,fp);
	    }
	}
      /* write Vel[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      Vel[0] = haloB->mhalos[ihalo].VXc;
	      Vel[1] = haloB->mhalos[ihalo].VYc;
	      Vel[2] = haloB->mhalos[ihalo].VZc;
	      fwrite(Vel,sizeof(float),3,fp);
	    }
	}
      /* write VelDisp */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      VelDisp = haloB->mhalos[ihalo].sigV;
	      fwrite(&(VelDisp),sizeof(float),1,fp);
	    }
	}     
      /* write Vmax */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      Vmax = haloB->mhalos[ihalo].Vmax;
	      fwrite(&(Vmax),sizeof(float),1,fp);
	    }
	}    
      /* write Spin[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      Spin[0] = haloB->mhalos[ihalo].SpinX;
	      Spin[1] = haloB->mhalos[ihalo].SpinY;
	      Spin[2] = haloB->mhalos[ihalo].SpinZ;
	      fwrite(Spin,sizeof(float),3,fp);
	    }
	}    
      /* write MostBoundID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  if(haloB->mhalos[ihalo].used == 1)
	    {
	      MostBoundID = haloB->mhalos[ihalo].Particles[0].ID;
	      fwrite(&(MostBoundID),sizeof(long long),1,fp);
	    }
	}    
      fclose(fp);
    }
  else
    {
      printf("Cannot open file %s\nExiting...\n",filename);
      exit(1);
    }
}


void internalaux_read(clgal_aux_data_wrapper_t *aux_data, char* outputfolder)
{

  hid_t ihalo,whalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  int ih;
  
  /* return if read */
  if(aux_data->already_read == 1)
    return;
  
  sprintf(filename,"%s/%3.3f/mtaux_%d.dat_bin",outputfolder,aux_data->redshift,aux_data->domainid);
  printf("reading %s\n",filename);
  fp = fopen(filename, "rb");
  if(fp != NULL)
    {
      /* read  total halo number */
      fread(&(aux_data->nHalos),sizeof(hid_t),1,fp);
      aux_data->lgal_aux_halos = calloc(aux_data->nHalos,sizeof(clgal_aux_data_t));
      /* read nprogs */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].nprogs),sizeof(uint32_t),1,fp);
	  aux_data->lgal_aux_halos[ihalo].proglist = calloc(aux_data->lgal_aux_halos[ihalo].nprogs,sizeof(hid_t));
	}
      /* read proglist */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  for(whalo=0; whalo < aux_data->lgal_aux_halos[ihalo].nprogs; whalo++)
	    {
	      fread(&(aux_data->lgal_aux_halos[ihalo].proglist[whalo]),sizeof(hid_t),1,fp);
	    }
	}

      /* write globalRefID */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].globalRefID),sizeof(hid_t),1,fp);
	}
      /* read FirstFOF globalRefID */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].FirstFOF),sizeof(hid_t),1,fp);
	}  
      /* read NextFOF globalRefID */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].NextFOF),sizeof(hid_t),1,fp);
	}  
      /* read M200 */
      
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.M_Crit200),sizeof(float),1,fp);
	  /* printf("Just read in : %llu --- %d : M:%f\n",aux_data->lgal_aux_halos[ihalo].globalRefID,aux_data->lgal_aux_halos[ihalo].nprogs,aux_data->lgal_aux_halos[ihalo].lgal_halo_data.M_Crit200); */
	  /* for(whalo=0;whalo<aux_data->lgal_aux_halos[ihalo].nprogs;whalo++) */
	  /*   printf("%llu: prog  ----> %llu\n",whalo,aux_data->lgal_aux_halos[ihalo].proglist[whalo]); */
	  //printf("M200 = %f\n",aux_data->lgal_aux_halos[ihalo].lgal_halo_data.M_Crit200);
	}
      /* read Len */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.Len),sizeof(int),1,fp);
	}
      /* read Pos[3] */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.Pos,sizeof(float),3,fp);
	}
      /* write Vel[3] */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.Vel,sizeof(float),3,fp);
	}
      /* write VelDisp */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.VelDisp),sizeof(float),1,fp);
	}     
      /* write Vmax */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.Vmax),sizeof(float),1,fp);
	}    
      /* write Spin[3] */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.Spin,sizeof(float),3,fp);
	}    
      /* write MostBoundID */
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  fread(&(aux_data->lgal_aux_halos[ihalo].lgal_halo_data.MostBoundID),sizeof(long long),1,fp);
	}    
      fclose(fp);
      for(ihalo=0; ihalo < aux_data->nHalos; ihalo++)
	{
	  if(aux_data->lgal_aux_halos[ihalo].FirstFOF == NULLPOINT)
	    {
	      aux_data->lgal_aux_halos[ihalo].FirstFOF = aux_data->lgal_aux_halos[ihalo].globalRefID;
	    }
	  aux_data->lgal_aux_halos[ihalo].FirstProgenitor = NULLPOINT;
	  aux_data->lgal_aux_halos[ihalo].NextProgenitor = NULLPOINT;
	  aux_data->lgal_aux_halos[ihalo].Descendant = NULLPOINT;
	  aux_data->lgal_aux_halos[ihalo].TreeNr = -1;
	  aux_data->lgal_aux_halos[ihalo].hidTree = -1;
	  aux_data->lgal_aux_halos[ihalo].doneaux = 0;
	  aux_data->lgal_aux_halos[ihalo].donetree = 0;
	  aux_data->lgal_aux_halos[ihalo].redshift = aux_data->redshift;
	}
      qsort(aux_data->lgal_aux_halos,aux_data->nHalos, sizeof(clgal_aux_data_t),compare_clgal_aux_data_t_by_globalRefID);

      aux_data->already_read = 1;

      /* if(aux_data->snapid > 0) */
      /* 	{ */
      /* 	  for(ihalo=0; ihalo < aux_data->nHalos; ihalo++) */
      /* 	    { */
      /* 	      printf("read in : %llu --- %d : M:%f\n",aux_data->lgal_aux_halos[ihalo].globalRefID,aux_data->lgal_aux_halos[ihalo].nprogs,aux_data->lgal_aux_halos[ihalo].lgal_halo_data.M_Crit200); */
      /* 	      for(whalo=0;whalo<aux_data->lgal_aux_halos[ihalo].nprogs;whalo++) */
      /* 		printf("%llu: prog  ----> %llu\n",whalo,aux_data->lgal_aux_halos[ihalo].proglist[whalo]); */
      /* 	      //printf("read in : %llu ---> %llu:%llu/%d\n",aux_data->lgal_aux_halos[ihalo].globalRefID,aux_data->lgal_aux_halos[ihalo].proglist[whalo],whalo,aux_data->lgal_aux_halos[ihalo].nprogs); */

      /* 	    } */
      /* 	} */
    }
  else
    {
      printf("Cannot open file %s\nExiting...\n",filename);
      exit(1);
    }

}

int compare_clgal_aux_data_t_by_globalRefID(const void *v1, const void *v2)
{
    const clgal_aux_data_t *u1 = v1;
    const clgal_aux_data_t *u2 = v2;
    int ret;
    if(u1->globalRefID < u2->globalRefID)
      ret =  -1;
    else if(u1->globalRefID > u2->globalRefID)
      ret = 1;
    else if(u1->globalRefID == u2->globalRefID)
      ret = 0;
    return ret;
}

uint64_t search_clgal_aux_data_t_for_globalRefID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high,ipart;
  clgal_aux_data_t *pool = (clgal_aux_data_t *) Array;

  low = 0;
  high = n_array-1;

  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
      //      printf("%llu : %llu < %llu : %llu < %llu : %llu\n",low,pool[low].ID,middle,pool[middle].ID,high,pool[high].ID);
      if ( searchID == pool[ middle ].globalRefID)
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].globalRefID)
	{
	  high = middle - 1;
	}
      else
	{
	  low = middle + 1;
	}
    }
  return NULLPOINT;
}


