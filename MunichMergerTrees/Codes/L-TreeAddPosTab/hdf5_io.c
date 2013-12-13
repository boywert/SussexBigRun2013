#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
void read_hdf5_attribute(hid_t file_id, 
			 char *obj_name, char *attr_name,
			 hid_t memtype_id, void *buf)
{
  hid_t attr_id;
  hid_t parent_id;

  parent_id = H5Gopen1(file_id, obj_name);
  if(parent_id<0)
    {
      printf("Failed to open HDF5 group\n");
      printf("Object name: %s\n", obj_name);
      fflush(stdout);
      exit(779);
    }
  attr_id = H5Aopen_name(parent_id, attr_name);
  if(attr_id < 0)
    {
      printf("Failed to open HDF5 attribute\n");
      printf("Object name: %s\n", obj_name);
      printf("Attribute name: %s\n", attr_name);
      fflush(stdout);
      exit(779);
    }
  if(H5Aread(attr_id, memtype_id, buf) < 0)
    {
      printf("Failed to read HDF5 attribute\n");
      printf("Object name: %s\n", obj_name);
      printf("Attribute name: %s\n", attr_name);
      fflush(stdout);
      exit(780);
    }
  H5Aclose(attr_id);
  H5Gclose(parent_id);
}

void read_hdf5_header(hid_t file_id, struct io_header *head)
{
  read_hdf5_attribute(file_id, "/Header", "NumPart_ThisFile",    H5T_NATIVE_INT,      head->npart);
  read_hdf5_attribute(file_id, "/Header", "MassTable",           H5T_NATIVE_DOUBLE,      head->mass);
  read_hdf5_attribute(file_id, "/Header", "Time",                H5T_NATIVE_DOUBLE, &(head->time));
  read_hdf5_attribute(file_id, "/Header", "Redshift",            H5T_NATIVE_DOUBLE, &(head->redshift));
  read_hdf5_attribute(file_id, "/Header", "Flag_Sfr",            H5T_NATIVE_INT,    &(head->flag_sfr));
  read_hdf5_attribute(file_id, "/Header", "Flag_Feedback",       H5T_NATIVE_INT,    &(head->flag_feedback));
  read_hdf5_attribute(file_id, "/Header", "NumPart_Total",       H5T_NATIVE_UINT,     head->npartTotal);
  read_hdf5_attribute(file_id, "/Header", "Flag_Cooling",        H5T_NATIVE_INT,    &(head->flag_cooling));
  read_hdf5_attribute(file_id, "/Header", "NumFilesPerSnapshot", H5T_NATIVE_INT,    &(head->num_files));
  read_hdf5_attribute(file_id, "/Header", "BoxSize",             H5T_NATIVE_DOUBLE, &(head->BoxSize));
  read_hdf5_attribute(file_id, "/Header", "Omega0",              H5T_NATIVE_DOUBLE, &(head->Omega0));
  read_hdf5_attribute(file_id, "/Header", "OmegaLambda",         H5T_NATIVE_DOUBLE, &(head->OmegaLambda));
  read_hdf5_attribute(file_id, "/Header", "HubbleParam",         H5T_NATIVE_DOUBLE, &(head->HubbleParam));
  read_hdf5_attribute(file_id, "/Header", "Flag_StellarAge",     H5T_NATIVE_INT,    &(head->flag_stellarage));
  read_hdf5_attribute(file_id, "/Header", "Flag_Metals",         H5T_NATIVE_INT,    &(head->flag_metals));
}

void read_hdf5_dataset(hid_t file_id, char *dset_name, hid_t dtype_id, void *buf)
{
  hid_t dset_id;

  if((dset_id = H5Dopen1(file_id, dset_name)) < 0)
    {
      printf("Failed to open HDF5 dataset\n");
      printf("Dataset name: %s\n", dset_name);
      fflush(stdout);
      exit(781);
    }
  if(H5Dread(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0)
    {
      printf("Failed to read HDF5 dataset\n");
      printf("Dataset name: %s\n", dset_name);
      fflush(stdout);
      exit(782);
    }
  H5Dclose(dset_id);
}

/*
  Read the quantity named by dset_name for all particle types in the file
  (e.g. dset_name can be "Coordinates", "Velocities", "ParticleIDs" etc)
*/
void read_all_types(hid_t file_id, struct io_header *head, char *dset_name, hid_t dtype_id, void *buf)
{
  hid_t dset_id, dspace_id;
  int itype;
  char name[500];
  char *ptr;
  hssize_t npoints;

  ptr = (char *) buf;
  for(itype=0; itype<6; itype+=1)
    {
      if(head->npart[itype] > 0)
	{
	  sprintf(name, "PartType%i/%s", itype, dset_name);
	  if((dset_id = H5Dopen1(file_id, name)) < 0)
	    {
	      printf("Failed to open HDF5 dataset\n");
	      printf("Dataset name: %s\n", name);
	      fflush(stdout);
	      exit(783); 
	    }
	  if(H5Dread(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr) < 0)
	    {
	      printf("Failed to read HDF5 dataset\n");
	      printf("Dataset name: %s\n", name);
	      fflush(stdout);
	      exit(784); 
	    }
	  dspace_id = H5Dget_space(dset_id);
	  npoints = H5Sget_simple_extent_npoints(dspace_id); 
	  ptr += (H5Tget_size(dtype_id) * npoints);
	  H5Sclose(dspace_id);
	  H5Dclose(dset_id);
	}
    }
} 
#endif
