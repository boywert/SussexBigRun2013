#include "main.h"

int main(int argc,char **argv)

{   
  char   outfile[2048], **infile;
  char   prefix[2048], suffix[2048], suffix_ascii[2048];
  char halofile[2048],partfile[2048],folder[1024];
  int    i, slen, ftype, nfiles;
  init_memmgr();

  sprintf(folder,"/export/research/virgo/Boyd/testcurie");

  sussexbigrun_load_halo_catalogue_binary(folder,6.418,6*6*6);
}
  
