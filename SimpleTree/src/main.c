#include "main.h"

int main(int argc,char **argv)

{   
  char folder[1024];
  m_halo_wrapper_t *halocatA,*halocatB;
  char memmgr_buff[memmgr_max_str];
  init_memmgr();
  //sprintf(folder,"/ccc/cont005/home/ra1089/schneida/scratch/AHF_haloes/cubepm_130329_10_4000_244Mpc_ext2_test/results");
  sprintf(memmgr_buff,"Halo wrapper");
  sprintf(folder,"/mnt/lustre/scratch/cs390/testcurie");  
  //halocat = sussexbigrun_load_halo_catalogue_binary(folder,6.000,10*10*10);
  halocatA = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  halocatB = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  printf("finish alloc\n");
  printf("NULLPOINT  = %lu, MAXUSE = %lu\n",(uint64_t)-1,(uint64_t)-2 );
  halocatA[0] = sussexbigrun_load_halo_catalogue_binary(folder,6.354,6);
  halocatA[0] = sussexbigrun_load_halo_catalogue_binary(folder,6.418,6);
  memmgr_printdetails();
  free_m_halo_wrapper(halocatA);
  memmgr_printdetails();
  return 0;
}
  
