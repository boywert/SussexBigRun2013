#include "common_io.h"

int ReadInt(FILE *fptr,int32_t *n,int swap)
{
  /* unsigned char *cptr,tmp; */
  
  /* if(sizeof(int) != 4) */
  /*  { */
  /*   fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int)); */
  /*   exit(0); */
  /*  } */
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  /* if (swap) { */
  /*   cptr = (unsigned char *)n; */
  /*   tmp     = cptr[0]; */
  /*   cptr[0] = cptr[3]; */
  /*   cptr[3] = tmp; */
  /*   tmp     = cptr[1]; */
  /*   cptr[1] = cptr[2]; */
  /*   cptr[2] = tmp; */
  /* } */
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,uint32_t *n,int swap)
{
  /* unsigned char *cptr,tmp; */
  
  /* if(sizeof(int) != 4) */
  /*  { */
  /*   fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int)); */
  /*   exit(0); */
  /*  } */
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  /* if (swap) { */
  /*   cptr = (unsigned char *)n; */
  /*   tmp     = cptr[0]; */
  /*   cptr[0] = cptr[3]; */
  /*   cptr[3] = tmp; */
  /*   tmp     = cptr[1]; */
  /*   cptr[1] = cptr[2]; */
  /*   cptr[2] = tmp; */
  /* } */
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned long integer
 */
int ReadULong(FILE *fptr,uint64_t *n,int swap)
{
  /*   unsigned char *cptr,tmp; */
  
  /* if(sizeof(unsigned long) == 4) */
  /*  { */
  /*   if (fread(n,4,1,fptr) != 1) */
  /*     return(FALSE); */
  /*   /\* if (swap) { *\/ */
  /*   /\*   cptr = (unsigned char *)n; *\/ */
  /*   /\*   tmp     = cptr[0]; *\/ */
  /*   /\*   cptr[0] = cptr[3]; *\/ */
  /*   /\*   cptr[3] = tmp; *\/ */
  /*   /\*   tmp     = cptr[1]; *\/ */
  /*   /\*   cptr[1] = cptr[2]; *\/ */
  /*   /\*   cptr[2] = tmp; *\/ */
  /*   /\* } *\/ */
  /*  } */
  /* else if(sizeof(unsigned long) == 8) */
  /*  { */
    if (fread(n,8,1,fptr) != 1)
      return(FALSE);
    /* if (swap) { */
    /*   cptr = (unsigned char *)n; */
    /*   tmp     = cptr[0]; */
    /*   cptr[0] = cptr[7]; */
    /*   cptr[7] = tmp; */
    /*   tmp     = cptr[1]; */
    /*   cptr[1] = cptr[6]; */
    /*   cptr[6] = tmp; */
    /*   tmp     = cptr[2]; */
    /*   cptr[2] = cptr[5]; */
    /*   cptr[5] = tmp; */
    /*   tmp     = cptr[3]; */
    /*   cptr[3] = cptr[4]; */
    /*   cptr[4] = tmp; */
    /* } */
  /*  } */
  /* else */
  /*  { */
  /*   fprintf(stderr,"ReadULong: something wrong...cannot read long\n"); */
  /*   exit(0); */
  /*  } */
  
  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr, float *n, int swap)
{
  /*  unsigned char *cptr,tmp; */
  
  /* if(sizeof(float) != 4) */
  /*  { */
  /*   fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float)); */
  /*   exit(0); */
  /*  } */
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  /* if (swap) */
  /*  { */
  /*   cptr = (unsigned char *)n; */
  /*   tmp     = cptr[0]; */
  /*   cptr[0] = cptr[3]; */
  /*   cptr[3] = tmp; */
  /*   tmp     = cptr[1]; */
  /*   cptr[1] = cptr[2]; */
  /*   cptr[2] = tmp; */
  /*  } */
  return(TRUE);
}

