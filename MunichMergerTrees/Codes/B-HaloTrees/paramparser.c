#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "paramparser.h"

typedef struct
  {
  int nparams;
  int reserved;
  char **key, **value;
  } params;

static void strtrim (char *string)
  {
  size_t len=strlen(string);
  size_t lw=0;
  while ((len>0) && ((string[len-1]==' ')
                  || (string[len-1]=='\t')
                  || (string[len-1]=='\n')))
    --len;
  string[len]=0;
  while ((string[lw]==' ') || (string[lw]=='\t') || (string[lw]=='\n'))
    ++lw;
  if (lw>0)
    memmove (string,string+lw,len+1-lw); 
  }

static void push_params (params *par, const char *key, const char *value)
  {
  if (par->reserved<=par->nparams)
    {
    par->reserved*=2;
    if (par->reserved == 0) par->reserved=1;
    par->key = realloc(par->key,(par->reserved)*sizeof(char *)); 
    par->value = realloc(par->value,(par->reserved)*sizeof(char *));
    }

  par->key[par->nparams] = malloc(strlen(key)+1);
  strcpy (par->key[par->nparams],key);
  par->value[par->nparams] = malloc(strlen(value)+1);
  strcpy (par->value[par->nparams],value);
  ++par->nparams;
  }

void *openParams (const char *file)
  {
  FILE *fd;
  params *par;
  int cnt=0;
  char line[2048];

  fd = fopen(file,"r");
  if (!fd) return NULL;

  par = malloc (sizeof(params));
  par->nparams=0;
  par->reserved=0;
  par->key=NULL;
  par->value=NULL;

  while (fgets (line, 2000, fd)!=NULL)
    {
    size_t len;
    strtrim(line);
    ++cnt;
    len=strlen(line);

    if (len>0)
      {
	//modification
	/*DELETEME  additional requirements:
	  allowing space as separator without '=': first separated string is key, all the rest is value
	*/
	if (line[0]!='#' && line[0] != '%')
	  {
	    char key[2048], value[2048];
	    size_t eqpos = 0;
	    
	    while ((eqpos<len) && (line[eqpos]!='=')) ++eqpos;
	    
	    if(eqpos==0)
	      {
		fprintf(stderr, "ERROR: no key specified at line %d:\n%s\n",cnt,line);
		fclose(fd); return NULL;
	      }
	    
	    if (eqpos==len)
	      {
		/*		fprintf(stderr, "WARNING: no \'=\' for separating key and value at line %d\n", cnt);	*/	
		eqpos = 0;
		while ((eqpos<len) && line[eqpos] != ' ' && line[eqpos] != '\t') ++eqpos;
	      }

	strncpy (key,line,eqpos);
	key[eqpos]=0;
	strcpy (value,line+eqpos+1);
	strtrim(key);
	strtrim(value);
	push_params(par,key,value);
      }
    }
  }
fclose(fd);
return par;
}

void closeParams (void *vpar)
  {
  params *par=vpar;
  if (par!=NULL)
    {
    int m;
    for (m=0; m<par->nparams; ++m)
      {
      free (par->key[m]);
      free (par->value[m]);
      }
    if (par->key!=NULL) { free (par->key); par->key=NULL; }
    if (par->value!=NULL) { free (par->value); par->value=NULL; }
    free (par);
    }
  }

int locateParam (void *vpar, const char *key)
  {
  params *par=vpar;
  int m;
  for (m=0; m<par->nparams; ++m)
    if (strcmp(key,par->key[m])==0) return m;
  return -1;
  }

int isVar(void* vpar,const char *key)
{
  if(locateParam(vpar, key) != -1)
    return 1; /* key found */
  else
    return 0; /* key not found */
}

int getVarInt (void *vpar, const char *key, int *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtol(par->value[idx],&endptr,10);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarLongLong (void *vpar, const char *key, long long *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtoll(par->value[idx],&endptr,10);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarDouble (void *vpar, const char *key, double *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtod(par->value[idx],&endptr);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarString(void* vpar, const char* key, char** result)
  {
  return getVarStringMem(vpar, key, result);
  }

int getVarStringMem (void *vpar, const char *key, char **result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  *result = malloc (strlen(par->value[idx])+1);
  if (!result) return 1;                   /* out of memory */
  strcpy (*result,par->value[idx]);
  return 0;
  }

int getVarStringProt (void *vpar, const char *key, char *result, int capacity)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  if (strlen(par->value[idx])>capacity-1) return 1;   /* result does not fit */
  strcpy (result,par->value[idx]);
  return 0;
  }

int getVarStringUnProt (void *vpar, const char *key, char *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  strcpy (result,par->value[idx]);
  return 0;
  }
