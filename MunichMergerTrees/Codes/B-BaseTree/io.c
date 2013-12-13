#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "allvars.h"
#include "proto.h"

/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
	  fflush(stdout);
	  exit(777);
	}
    }
  else
    nwritten = 0;

  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) has occured: end of file\n");
      else
	printf("I/O error (fread) has occured: %s\n", strerror(errno));
      fflush(stdout);
      exit(778);
    }
  return nread;
}


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
  read_hdf5_attribute(file_id, "/Header", "NumPart_Total",       H5T_NATIVE_UINT,     head->npartTotal);
#ifdef HAVE_NUMPART_HIGHWORD
  read_hdf5_attribute(file_id, "/Header", "NumPart_Total_HighWord",H5T_NATIVE_UINT,   head->npartTotalHighWord);
#else
  int i;
  for(i=0;i<6;i+=1)
    head->npartTotalHighWord[i] = 0;
#endif
  read_hdf5_attribute(file_id, "/Header", "MassTable",           H5T_NATIVE_DOUBLE,      head->mass);
  read_hdf5_attribute(file_id, "/Header", "Time",                H5T_NATIVE_DOUBLE, &(head->time));
  read_hdf5_attribute(file_id, "/Header", "Redshift",            H5T_NATIVE_DOUBLE, &(head->redshift));
  read_hdf5_attribute(file_id, "/Header", "BoxSize",             H5T_NATIVE_DOUBLE, &(head->BoxSize));
  read_hdf5_attribute(file_id, "/Header", "NumFilesPerSnapshot", H5T_NATIVE_INT,    &(head->num_files));
  read_hdf5_attribute(file_id, "/Header", "Omega0",              H5T_NATIVE_DOUBLE, &(head->Omega0));
  read_hdf5_attribute(file_id, "/Header", "OmegaLambda",         H5T_NATIVE_DOUBLE, &(head->OmegaLambda));
  read_hdf5_attribute(file_id, "/Header", "HubbleParam",         H5T_NATIVE_DOUBLE, &(head->HubbleParam));
  read_hdf5_attribute(file_id, "/Header", "Flag_Sfr",            H5T_NATIVE_INT,    &(head->flag_sfr));
  read_hdf5_attribute(file_id, "/Header", "Flag_Cooling",        H5T_NATIVE_INT,    &(head->flag_cooling));
  read_hdf5_attribute(file_id, "/Header", "Flag_StellarAge",     H5T_NATIVE_INT,    &(head->flag_stellarage));
  read_hdf5_attribute(file_id, "/Header", "Flag_Metals",         H5T_NATIVE_INT,    &(head->flag_metals));
  read_hdf5_attribute(file_id, "/Header", "Flag_Feedback",       H5T_NATIVE_INT,    &(head->flag_feedback));
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
#endif
