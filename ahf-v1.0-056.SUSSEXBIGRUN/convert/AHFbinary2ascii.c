#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>

//=================================================================================
// This code aims at converting (and merging) binary _halos, _profiles,
// and _particles files generated by AHF when switching on -DAHFbinary
//
// This code also gives you the option to only convert those halo properties
// you are primarily interested in. But for that you need to adjust the
// part where the information is written commenting out those properties
// of no concern to you. Please adjust the two routines
//                   write_halos()   &   write_profiles()
// to best suit your needs. You are also able to filter out halos below a certain
// MassCut given on the commandline when starting this code.
//
// Note, none of the filtering works for _particles as this file does not contain
// any information about halo properties! Or to be more presice, fiddling with the
// _particles files is presently not supported.
//
// THIS CODE IS INDEPENDENT OF ANY AHF LIBRARY AND COULD BE COMPILED STAND-ALONE
//=================================================================================

#define VERBOSE
#define PROFILE_FOR_GNUPLOT  // will have blank lines separating each profile

//===================================================================
// the headers for the ASCII output files
//===================================================================
#define HEADER_HALOS    "# up to you what to put here\n"
#define HEADER_PROFILES "# up to you what to put here\n"

//===================================================================
// global variables
//===================================================================
double MassCut; // haloes below this mass will not be written to ASCII

//===================================================================
// structures
//===================================================================
typedef struct halo {
  uint32_t numColumns;
  uint64_t ID;
  uint64_t hostHalo;
  uint32_t numSubStruct;
  float    Mvir;
  uint32_t npart;
  float    Xc;
  float    Yc;
  float    Zc;
  float    VXc;
  float    VYc;
  float    VZc;
  float    Rvir;
  float    Rmax;
  float    r2;
  float    mbp_offset;
  float    com_offset;
  float    Vmax;
  float    v_esc;
  float    sigV;
  float    lambda;
  float    lambdaE;
  float    Lx;
  float    Ly;
  float    Lz;
  float    b;
  float    c;
  float    Eax;
  float    Eay;
  float    Eaz;
  float    Ebx;
  float    Eby;
  float    Ebz;
  float    Ecx;
  float    Ecy;
  float    Ecz;
  float    ovdens;
  uint32_t nbins;
  float    fMhires;
  float    Ekin;
  float    Epot;
  float    SurfP;
  float    Phi0;
  float    cNFW;
  uint32_t n_gas;
  float    M_gas;
  float    lambda_gas;
  float    lambdaE_gas;
  float    Lx_gas;
  float    Ly_gas;
  float    Lz_gas;
  float    b_gas;
  float    c_gas;
  float    Eax_gas;
  float    Eay_gas;
  float    Eaz_gas;
  float    Ebx_gas;
  float    Eby_gas;
  float    Ebz_gas;
  float    Ecx_gas;
  float    Ecy_gas;
  float    Ecz_gas;
  float    Ekin_gas;
  float    Epot_gas;
  uint32_t n_star;
  float    M_star;
  float    lambda_star;
  float    lambdaE_star;
  float    Lx_star;
  float    Ly_star;
  float    Lz_star;
  float    b_star;
  float    c_star;
  float    Eax_star;
  float    Eay_star;
  float    Eaz_star;
  float    Ebx_star;
  float    Eby_star;
  float    Ebz_star;
  float    Ecx_star;
  float    Ecy_star;
  float    Ecz_star;
  float    Ekin_star;
  float    Epot_star;
  float    mean_z_gas;
  float    mean_z_star;
} halo_t;

typedef struct halo_profile {
  uint32_t      nbins;
  uint32_t      numColumns;
  float         *r;
  uint32_t      *npart;
  float         *M_in_r;
  float         *ovdens;
  float         *dens;
  float         *vcirc;
  float         *vesc;
  float         *sigv;
  float         *Lx;
  float         *Ly;
  float         *Lz;
  float         *b;
  float         *c;
  float         *Eax;
  float         *Eay;
  float         *Eaz;
  float         *Ebx;
  float         *Eby;
  float         *Ebz;
  float         *Ecx;
  float         *Ecy;
  float         *Ecz;
  float         *Ekin;
  float         *Epot;
  float         *M_gas;
  float         *M_star;
  float         *u_gas;
  float         *Z_gas_sh;
  float         *Z_star_sh;
} halo_profile_t;

//===================================================================
// PROTOTYPES
//===================================================================
void     convert                           (int , int , char **, char *);
void     convert_halos                     (FILE *, FILE *);
void     convert_profiles                  (FILE *, FILE *);
void     convert_particles                 (FILE *, FILE *);
void     alloc_profile                     (uint64_t, halo_profile_t *);
void     free_profile                      (halo_profile_t *);
uint64_t get_TotalNumHalos_from_particles  (int, char**);

//===================================================================
// routines copied from AHF to keep this code independent
//===================================================================
int     ReadULong          (FILE *, unsigned long *, int);
int     ReadInt            (FILE *, int           *, int);
int     ReadUInt           (FILE *, unsigned int  *, int);
int     ReadFloat          (FILE *, float         *, int);
#define TRUE  1
#define FALSE 0

/*=================================================================================================
 * write_halo(): FEEL FREE TO FILTER UNWANTED PROPERTIES BY COMMENTING THEM OUT
 *=================================================================================================*/
void write_halo(halo_t halo, FILE *fpout)
{
  if(halo.Mvir >= MassCut){
    fprintf(fpout,"%ld ",(long unsigned)halo.ID);
    fprintf(fpout,"%ld ",(long unsigned)halo.hostHalo);
    fprintf(fpout,"%d ",halo.numSubStruct);
    fprintf(fpout,"%g ",halo.Mvir);
    fprintf(fpout,"%d ",halo.npart);
    fprintf(fpout,"%f ",halo.Xc);
    fprintf(fpout,"%f ",halo.Yc);
    fprintf(fpout,"%f ",halo.Zc);
    fprintf(fpout,"%f ",halo.VXc);
    fprintf(fpout,"%f ",halo.VYc);
    fprintf(fpout,"%f ",halo.VZc);
    fprintf(fpout,"%f ",halo.Rvir);
    fprintf(fpout,"%f ",halo.Rmax);
    fprintf(fpout,"%f ",halo.r2);
    fprintf(fpout,"%f ",halo.mbp_offset);
    fprintf(fpout,"%f ",halo.com_offset);
    fprintf(fpout,"%f ",halo.Vmax);
    fprintf(fpout,"%f ",halo.v_esc);
    fprintf(fpout,"%f ",halo.sigV);
    fprintf(fpout,"%f ",halo.lambda);
    fprintf(fpout,"%f ",halo.lambdaE);
    fprintf(fpout,"%f ",halo.Lx);
    fprintf(fpout,"%f ",halo.Ly);
    fprintf(fpout,"%f ",halo.Lz);
    fprintf(fpout,"%f ",halo.b);
    fprintf(fpout,"%f ",halo.c);
    fprintf(fpout,"%f ",halo.Eax);
    fprintf(fpout,"%f ",halo.Eay);
    fprintf(fpout,"%f ",halo.Eaz);
    fprintf(fpout,"%f ",halo.Ebx);
    fprintf(fpout,"%f ",halo.Eby);
    fprintf(fpout,"%f ",halo.Ebz);
    fprintf(fpout,"%f ",halo.Ecx);
    fprintf(fpout,"%f ",halo.Ecy);
    fprintf(fpout,"%f ",halo.Ecz);
    fprintf(fpout,"%f ",halo.ovdens);
    fprintf(fpout,"%d ",halo.nbins);
    fprintf(fpout,"%f ",halo.fMhires);
    fprintf(fpout,"%g ",halo.Ekin);
    fprintf(fpout,"%g ",halo.Epot);
    fprintf(fpout,"%g ",halo.SurfP);
    fprintf(fpout,"%g ",halo.Phi0);
    fprintf(fpout,"%f ",halo.cNFW);
    if(halo.numColumns > 43) {
      fprintf(fpout,"%d ",halo.n_gas);
      fprintf(fpout,"%g ",halo.M_gas);
      fprintf(fpout,"%f ",halo.lambda_gas);
      fprintf(fpout,"%f ",halo.lambdaE_gas);
      fprintf(fpout,"%f ",halo.Lx_gas);
      fprintf(fpout,"%f ",halo.Ly_gas);
      fprintf(fpout,"%f ",halo.Lz_gas);
      fprintf(fpout,"%f ",halo.b_gas);
      fprintf(fpout,"%f ",halo.c_gas);
      fprintf(fpout,"%f ",halo.Eax_gas);
      fprintf(fpout,"%f ",halo.Eay_gas);
      fprintf(fpout,"%f ",halo.Eaz_gas);
      fprintf(fpout,"%f ",halo.Ebx_gas);
      fprintf(fpout,"%f ",halo.Eby_gas);
      fprintf(fpout,"%f ",halo.Ebz_gas);
      fprintf(fpout,"%f ",halo.Ecx_gas);
      fprintf(fpout,"%f ",halo.Ecy_gas);
      fprintf(fpout,"%f ",halo.Ecz_gas);
      fprintf(fpout,"%g ",halo.Ekin_gas);
      fprintf(fpout,"%g ",halo.Epot_gas);
      fprintf(fpout,"%d ",halo.n_star);
      fprintf(fpout,"%g ",halo.M_star);
      fprintf(fpout,"%f ",halo.lambda_star);
      fprintf(fpout,"%f ",halo.lambdaE_star);
      fprintf(fpout,"%f ",halo.Lx_star);
      fprintf(fpout,"%f ",halo.Ly_star);
      fprintf(fpout,"%f ",halo.Lz_star);
      fprintf(fpout,"%f ",halo.b_star);
      fprintf(fpout,"%f ",halo.c_star);
      fprintf(fpout,"%f ",halo.Eax_star);
      fprintf(fpout,"%f ",halo.Eay_star);
      fprintf(fpout,"%f ",halo.Eaz_star);
      fprintf(fpout,"%f ",halo.Ebx_star);
      fprintf(fpout,"%f ",halo.Eby_star);
      fprintf(fpout,"%f ",halo.Ebz_star);
      fprintf(fpout,"%f ",halo.Ecx_star);
      fprintf(fpout,"%f ",halo.Ecy_star);
      fprintf(fpout,"%f ",halo.Ecz_star);
      fprintf(fpout,"%g ",halo.Ekin_star);
      fprintf(fpout,"%g ",halo.Epot_star);
    }
    if(halo.numColumns > 83) {
      fprintf(fpout,"%f ",halo.mean_z_gas);
      fprintf(fpout,"%f ",halo.mean_z_star);
    }
    fprintf(fpout,"\n");
  } // if(MassCut)
}

/*=================================================================================================
 * write_profile(): FEEL FREE TO FILTER UNWANTED PROPERTIES BY COMMENTING THEM OUT
 *=================================================================================================*/
void write_profile(halo_profile_t profile, FILE *fpout)
{
  uint32_t ibin;
  
  if(profile.M_in_r[profile.nbins-1] >= MassCut)
   {
#ifdef VERBOSE2
    fprintf(stderr,"   - halo id %ld (bins=%d)\n",(long unsigned)i,profile.nbins);
#endif
    for(ibin=0; ibin<profile.nbins; ibin++) {
      fprintf(fpout,"%f ", profile.r[ibin]);
      fprintf(fpout,"%d ", (int)profile.npart[ibin]);
      fprintf(fpout,"%g ", profile.M_in_r[ibin]);
      fprintf(fpout,"%f ", profile.ovdens[ibin]);
      fprintf(fpout,"%f ", profile.dens[ibin]);
      fprintf(fpout,"%f ", profile.vcirc[ibin]);
      fprintf(fpout,"%f ", profile.vesc[ibin]);
      fprintf(fpout,"%f ", profile.sigv[ibin]);
      fprintf(fpout,"%g ", profile.Lx[ibin]);
      fprintf(fpout,"%g ", profile.Ly[ibin]);
      fprintf(fpout,"%g ", profile.Lz[ibin]);
      fprintf(fpout,"%f ", profile.b[ibin]);
      fprintf(fpout,"%f ", profile.c[ibin]);
      fprintf(fpout,"%f ", profile.Eax[ibin]);
      fprintf(fpout,"%f ", profile.Eay[ibin]);
      fprintf(fpout,"%f ", profile.Eaz[ibin]);
      fprintf(fpout,"%f ", profile.Ebx[ibin]);
      fprintf(fpout,"%f ", profile.Eby[ibin]);
      fprintf(fpout,"%f ", profile.Ebz[ibin]);
      fprintf(fpout,"%f ", profile.Ecx[ibin]);
      fprintf(fpout,"%f ", profile.Ecy[ibin]);
      fprintf(fpout,"%f ", profile.Ecz[ibin]);
      fprintf(fpout,"%g ", profile.Ekin[ibin]);
      fprintf(fpout,"%g ", profile.Epot[ibin]);
      if(profile.numColumns > 24) {
        fprintf(fpout,"%g ", profile.M_gas[ibin]);
        fprintf(fpout,"%g ", profile.M_star[ibin]);
        fprintf(fpout,"%f ", profile.u_gas[ibin]);
      }
      if(profile.numColumns > 27) {
        fprintf(fpout,"%f ", profile.Z_gas_sh[ibin]);
        fprintf(fpout,"%f ", profile.Z_star_sh[ibin]);
      }
      fprintf(fpout,"\n");
    }
#ifdef PROFILE_FOR_GNUPLOT
    fprintf(fpout,"\n");
#endif
   } // if(MassCut)
}

/*=================================================================================================
 * main()
 *=================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{   
  char   outfile[2048], **infile;
  char   prefix[2048], suffix[2048], suffix_ascii[2048];

  int    i, slen, ftype, nfiles;

  fprintf(stderr,"======================================================================\n");
  fprintf(stderr," convert files generated by AHF when using -DAHFbinary to ASCII files\n");
  fprintf(stderr,"   (the code also concatenates all catalogues when run in MPI mode)\n");
  fprintf(stderr,"======================================================================\n");

  if(argc<5)
    {
      fprintf(stderr,"usage: %s prefix suffix Nfiles MassCut\n",*argv);
      exit(1);
    }
  
  //===================================================================
  // deal with command line
  //===================================================================
  strcpy(prefix,argv[1]);
  strcpy(suffix,argv[2]);
  nfiles  = atoi(argv[3]);
  MassCut = atof(argv[4]);
  
  // figure out file type
  ftype = -1;
  slen  = strlen(suffix);
  if(strcmp(suffix+(slen-9), "halos_bin") == 0) {
    ftype = 0;
    fprintf(stderr,"o will be converting %d halos files\n",nfiles);
  }
  if(strcmp(suffix+(slen-12), "profiles_bin") == 0) {
    ftype = 1;
    fprintf(stderr,"o will be converting %d profiles files\n",nfiles);
  }
  if(strcmp(suffix+(slen-13), "particles_bin") == 0) {
    ftype = 2;
    fprintf(stderr,"o will be converting %d particles files\n",nfiles);
  }

  // prepare array holding filenames
  infile = (char **) calloc(nfiles, sizeof(char *));
  for (i=0; i<nfiles; i++) {
    infile[i] = (char *) calloc(2048, sizeof(char));
  }

  // construct filenames
  if(nfiles == 1) {
    sprintf(infile[0],"%s.%s",prefix,suffix);
  }
  else {
    for(i=0; i<nfiles; i++) {
      sprintf(infile[i],"%s.%04d.%s",prefix,i,suffix);
    }
  }
  strncpy(suffix_ascii,suffix,slen-4);
  //suffix_ascii[slen-3] = "\0";
  sprintf(outfile,"%s.%s",prefix,suffix_ascii);
  
  // be verbose
  fprintf(stderr,"o freading from:\n");
  for(i=0; i<nfiles; i++) {
    fprintf(stderr,"   %s\n",infile[i]);
  }
  fprintf(stderr,"o writing to: %s\n",outfile);
  
  //===================================================================
  // wrapper for file-type conversion routines
  //===================================================================
  convert(ftype, nfiles, infile, outfile);
}
  
/*=================================================================================================
 * convert_halos()
 *=================================================================================================*/
void  convert(int ftype, int nfiles, char **infile, char *outfile)
{
  int i;
  FILE *fpin, *fpout;
  uint64_t TotalNumHalos;
  
  // open output file
  fpout = fopen(outfile,"w");
  if(fpout == NULL) {
    fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
    exit(0);
  }
  
  // write the header line
  switch (ftype) {
    case 0:
      fprintf(fpout,HEADER_HALOS);
      break;
    case 1:
      fprintf(fpout,HEADER_PROFILES);
      break;
    case 2:
      TotalNumHalos = get_TotalNumHalos_from_particles(nfiles,infile);
      fprintf(fpout,"%ld\n",(long unsigned)TotalNumHalos);
      break;
  }
  
  // loop over all input files
  for(i=0; i<nfiles; i++) {
    
    // open input file
    fpin = fopen(infile[i],"rb");
    if(fpin == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for writing\n",infile[i]);
      exit(0);
    }
    
    // which format to convert
    switch (ftype) {
      case 0:
        convert_halos(fpin,fpout);
        break;
      case 1:
        convert_profiles(fpin,fpout);
        break;
      case 2:
        convert_particles(fpin,fpout);
        break;
    }
    fclose(fpin);
  }
  
  fclose(fpout);
}

/*=================================================================================================
 * convert_halos()
 *=================================================================================================*/
void  convert_halos(FILE *fpin, FILE *fpout)
{
  uint64_t numHalos;
  uint32_t numColumns;
  uint64_t i;
  int32_t  one;
  int      swap=0;
  halo_t   halo;

  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)    swap = 0;
  else            swap = 1;

  ReadULong(fpin, &numHalos,   swap);
  ReadUInt (fpin, &numColumns, swap);
  
#ifdef VERBOSE
  fprintf(stderr,"o reading %ld halos from file (swap=%d,numColumns=%d)\n",(long unsigned)numHalos,swap,numColumns);
#endif

  // read in halo properties
  for(i=0; i<numHalos; i++) {
    ReadULong(fpin, &halo.ID,           swap);    // ID(1)
    ReadULong(fpin, &halo.hostHalo,     swap);    // hostHalo(2)
    ReadUInt (fpin, &halo.numSubStruct, swap);    // numSubStruct(3)
    ReadFloat(fpin, &halo.Mvir,         swap);    // Mvir(4)
    ReadUInt (fpin, &halo.npart,        swap);    // npart(5)
    ReadFloat(fpin, &halo.Xc,           swap);    // Xc(6)
    ReadFloat(fpin, &halo.Yc,           swap);    // Yc(7)
    ReadFloat(fpin, &halo.Zc,           swap);    // Zc(8)
    ReadFloat(fpin, &halo.VXc,          swap);    // VXc(9)
    ReadFloat(fpin, &halo.VYc,          swap);    // VYc(10)
    ReadFloat(fpin, &halo.VZc,          swap);    // VZc(11)
    ReadFloat(fpin, &halo.Rvir,         swap);    // Rvir(12)
    ReadFloat(fpin, &halo.Rmax,         swap);    // Rmax(13)
    ReadFloat(fpin, &halo.r2,           swap);    // r2(14)
    ReadFloat(fpin, &halo.mbp_offset,   swap);    // mbp_offset(15)
    ReadFloat(fpin, &halo.com_offset,   swap);    // com_offset(16)
    ReadFloat(fpin, &halo.Vmax,         swap);    // Vmax(17)
    ReadFloat(fpin, &halo.v_esc,        swap);    // v_esc(18)
    ReadFloat(fpin, &halo.sigV,         swap);    // sigV(19)
    ReadFloat(fpin, &halo.lambda,       swap);    // lambda(20)
    ReadFloat(fpin, &halo.lambdaE,      swap);    // lambdaE(21)
    ReadFloat(fpin, &halo.Lx,           swap);    // Lx(22)
    ReadFloat(fpin, &halo.Ly,           swap);    // Ly(23)
    ReadFloat(fpin, &halo.Lz,           swap);    // Lz(24)
    ReadFloat(fpin, &halo.b,            swap);    // b(25)
    ReadFloat(fpin, &halo.c,            swap);    // c(26)
    ReadFloat(fpin, &halo.Eax,          swap);    // Eax(27)
    ReadFloat(fpin, &halo.Eay,          swap);    // Eay(28)
    ReadFloat(fpin, &halo.Eaz,          swap);    // Eaz(29)
    ReadFloat(fpin, &halo.Ebx,          swap);    // Ebx(30)
    ReadFloat(fpin, &halo.Eby,          swap);    // Eby(31)
    ReadFloat(fpin, &halo.Ebz,          swap);    // Ebz(32)
    ReadFloat(fpin, &halo.Ecx,          swap);    // Ecx(33)
    ReadFloat(fpin, &halo.Ecy,          swap);    // Ecy(34)
    ReadFloat(fpin, &halo.Ecz,          swap);    // Ecz(35)
    ReadFloat(fpin, &halo.ovdens,       swap);    // ovdens(36)
    ReadUInt (fpin, &halo.nbins,        swap);    // nbins(37)
    ReadFloat(fpin, &halo.fMhires,      swap);    // fMhires(38)
    ReadFloat(fpin, &halo.Ekin,         swap);    // Ekin(39)
    ReadFloat(fpin, &halo.Epot,         swap);    // Epot(40)
    ReadFloat(fpin, &halo.SurfP,        swap);    // SurfP(41)
    ReadFloat(fpin, &halo.Phi0,         swap);    // Phi0(42)
    ReadFloat(fpin, &halo.cNFW,         swap);    // cNFW(43)
    if(numColumns > 43) {
      ReadUInt (fpin, &halo.n_gas,       swap);    // n_gas(44)
      ReadFloat(fpin, &halo.M_gas,       swap);    // M_gas(45)
      ReadFloat(fpin, &halo.lambda_gas,  swap);    // lambda_gas(46)
      ReadFloat(fpin, &halo.lambdaE_gas, swap);    // lambdaE_gas(47)
      ReadFloat(fpin, &halo.Lx_gas,      swap);    // Lx_gas(48)
      ReadFloat(fpin, &halo.Ly_gas,      swap);    // Ly_gas(49)
      ReadFloat(fpin, &halo.Lz_gas,      swap);    // Lz_gas(50)
      ReadFloat(fpin, &halo.b_gas,       swap);    // b_gas(51)
      ReadFloat(fpin, &halo.c_gas,       swap);    // c_gas(52)
      ReadFloat(fpin, &halo.Eax_gas,     swap);    // Eax_gas(53)
      ReadFloat(fpin, &halo.Eay_gas,     swap);    // Eay_gas(54)
      ReadFloat(fpin, &halo.Eaz_gas,     swap);    // Eaz_gas(55)
      ReadFloat(fpin, &halo.Ebx_gas,     swap);    // Ebx_gas(56)
      ReadFloat(fpin, &halo.Eby_gas,     swap);    // Eby_gas(57)
      ReadFloat(fpin, &halo.Ebz_gas,     swap);    // Ebz_gas(58)
      ReadFloat(fpin, &halo.Ecx_gas,     swap);    // Ecx_gas(59)
      ReadFloat(fpin, &halo.Ecy_gas,     swap);    // Ecy_gas(60)
      ReadFloat(fpin, &halo.Ecz_gas,     swap);    // Ecz_gas(61)
      ReadFloat(fpin, &halo.Ekin_gas,    swap);    // Ekin_gas(62)
      ReadFloat(fpin, &halo.Epot_gas,    swap);    // Epot_gas(63)
      ReadUInt (fpin, &halo.n_star,       swap);    // n_star(64)
      ReadFloat(fpin, &halo.M_star,       swap);    // M_star(65)
      ReadFloat(fpin, &halo.lambda_star,  swap);    // lambda_star(66)
      ReadFloat(fpin, &halo.lambdaE_star, swap);    // lambdaE_star(67)
      ReadFloat(fpin, &halo.Lx_star,      swap);    // Lx_star(68)
      ReadFloat(fpin, &halo.Ly_star,      swap);    // Ly_star(69)
      ReadFloat(fpin, &halo.Lz_star,      swap);    // Lz_star(70)
      ReadFloat(fpin, &halo.b_star,       swap);    // b_star(71)
      ReadFloat(fpin, &halo.c_star,       swap);    // c_star(72)
      ReadFloat(fpin, &halo.Eax_star,     swap);    // Eax_star(73)
      ReadFloat(fpin, &halo.Eay_star,     swap);    // Eay_star(74)
      ReadFloat(fpin, &halo.Eaz_star,     swap);    // Eaz_star(75)
      ReadFloat(fpin, &halo.Ebx_star,     swap);    // Ebx_star(76)
      ReadFloat(fpin, &halo.Eby_star,     swap);    // Eby_star(77)
      ReadFloat(fpin, &halo.Ebz_star,     swap);    // Ebz_star(78)
      ReadFloat(fpin, &halo.Ecx_star,     swap);    // Ecx_star(79)
      ReadFloat(fpin, &halo.Ecy_star,     swap);    // Ecy_star(80)
      ReadFloat(fpin, &halo.Ecz_star,     swap);    // Ecz_star(81)
      ReadFloat(fpin, &halo.Ekin_star,    swap);    // Ekin_star(82)
      ReadFloat(fpin, &halo.Epot_star,    swap);    // Epot_star(83)
    }
    if(numColumns > 83) {
      ReadFloat(fpin, &halo.mean_z_gas,    swap);    // mean_z_gas(84)
      ReadFloat(fpin, &halo.mean_z_star,   swap);    // mean_z_star(85)
    }
    
    //=================================================================================
    // write halo properties
    //=================================================================================
    halo.numColumns = numColumns;
    write_halo(halo, fpout);
    
  } // for(numHalos)
}

/*=================================================================================================
 * convert_profiles()
 *=================================================================================================*/
void  convert_profiles(FILE *fpin, FILE *fpout)
{
  uint64_t       numHalos;
  uint32_t       numColumns;
  uint64_t       i;
  uint32_t       ibin, nbins;
  halo_profile_t profile;
  int32_t        one;
  int            swap;

  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)    swap = 0;
  else            swap = 1;
  
  ReadULong(fpin, &numHalos,   swap);
  ReadUInt (fpin, &numColumns, swap);

#ifdef VERBOSE
  fprintf(stderr,"o reading %ld halos from file (swap=%d,numColumns=%d)\n",(long unsigned)numHalos,swap,numColumns);
#endif

  for(i=0; i<numHalos; i++) {
    
    // read number of bins
    ReadUInt(fpin, &nbins, swap);
    
#ifdef VERBOSE2
    fprintf(stderr,"   - halo id %ld (bins=%d)\n",(long unsigned)i,nbins);
#endif
    
    // allocate profile bins
    alloc_profile(nbins, &profile);
    
    // read full profile in
    for(ibin=0; ibin<nbins; ibin++) {
      ReadFloat(fpin, &(profile.r[ibin]),      swap);
      ReadUInt (fpin, &(profile.npart[ibin]),  swap);
      ReadFloat(fpin, &(profile.M_in_r[ibin]), swap);
      ReadFloat(fpin, &(profile.ovdens[ibin]), swap);
      ReadFloat(fpin, &(profile.dens[ibin]),   swap);
      ReadFloat(fpin, &(profile.vcirc[ibin]),  swap);
      ReadFloat(fpin, &(profile.vesc[ibin]),   swap);
      ReadFloat(fpin, &(profile.sigv[ibin]),   swap);
      ReadFloat(fpin, &(profile.Lx[ibin]),     swap);
      ReadFloat(fpin, &(profile.Ly[ibin]),     swap);
      ReadFloat(fpin, &(profile.Lz[ibin]),     swap);
      ReadFloat(fpin, &(profile.b[ibin]),      swap);
      ReadFloat(fpin, &(profile.c[ibin]),      swap);
      ReadFloat(fpin, &(profile.Eax[ibin]),    swap);
      ReadFloat(fpin, &(profile.Eay[ibin]),    swap);
      ReadFloat(fpin, &(profile.Eaz[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ebx[ibin]),    swap);
      ReadFloat(fpin, &(profile.Eby[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ebz[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ecx[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ecy[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ecz[ibin]),    swap);
      ReadFloat(fpin, &(profile.Ekin[ibin]),   swap);
      ReadFloat(fpin, &(profile.Epot[ibin]),   swap);
      if(numColumns > 24) {
        ReadFloat(fpin, &(profile.M_gas[ibin]),   swap);
        ReadFloat(fpin, &(profile.M_star[ibin]),  swap);
        ReadFloat(fpin, &(profile.u_gas[ibin]),   swap);
      }
      if(numColumns > 27) {
        ReadFloat(fpin, &(profile.Z_gas_sh[ibin]),   swap);
        ReadFloat(fpin, &(profile.Z_star_sh[ibin]),  swap);
      }
    }
    
    //=================================================================================
    // write halo profile
    //=================================================================================
    profile.nbins      = nbins;
    profile.numColumns = numColumns;
    write_profile(profile, fpout);
    
    // remove profile bins again
    free_profile(&profile);
    
  } // for(numHalos)
  

}

/*=================================================================================================
 * convert_particles()
 *=================================================================================================*/
void  convert_particles(FILE *fpin, FILE *fpout)
{
  uint64_t       numHalos, ihalo, npart, ipart, id;
  uint32_t       numColumns;
  int32_t        ptype;
  int32_t        one;
  int            swap;
#ifdef SUSSEXBIGRUN
  float          energy;
#endif
  
  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)   swap = 0;
  else           swap = 1;
  
  ReadULong(fpin, &numHalos,   swap);
  ReadUInt (fpin, &numColumns, swap);
  
#ifdef VERBOSE
  fprintf(stderr,"o reading %ld halos from file (swap=%d)\n",(long unsigned)numHalos,swap);
#endif

  // as we are not filtering anything we read & write as we go
  for(ihalo=0; ihalo<numHalos; ihalo++) {
    ReadULong(fpin, &npart, swap);
    fprintf(fpout,"%ld\n",(long unsigned)npart);
    for(ipart=0; ipart<npart; ipart++) {
      ReadULong(fpin, &id, swap);
      fprintf(fpout,"%ld ",(long unsigned)id);
      if(numColumns>1) {
#ifdef SUSSEXBIGRUN
	ReadFloat(fpin, &energy, swap);
	fprintf(fpout, "%g", (float)energy);
#else
        ReadInt(fpin, &ptype, swap);
        fprintf(fpout,"%d",(int)ptype);
#endif
      }
      fprintf(fpout,"\n");
    }
  }
}


/*=================================================================================================
 * alloc_profile()
 *=================================================================================================*/
void alloc_profile(uint64_t nbins, halo_profile_t *prof)
{
  prof->r       = (float *)       calloc(nbins, sizeof(float));
  prof->npart   = (uint32_t*)     calloc(nbins, sizeof(uint32_t));
  prof->M_in_r  = (float *)       calloc(nbins, sizeof(float));
  prof->ovdens  = (float *)       calloc(nbins, sizeof(float));
  prof->dens    = (float *)       calloc(nbins, sizeof(float));
  prof->vcirc   = (float *)       calloc(nbins, sizeof(float));
  prof->vesc    = (float *)       calloc(nbins, sizeof(float));
  prof->sigv    = (float *)       calloc(nbins, sizeof(float));
  prof->Lx      = (float *)       calloc(nbins, sizeof(float));
  prof->Ly      = (float *)       calloc(nbins, sizeof(float));
  prof->Lz      = (float *)       calloc(nbins, sizeof(float));
  prof->b       = (float *)       calloc(nbins, sizeof(float));
  prof->c       = (float *)       calloc(nbins, sizeof(float));
  prof->Eax     = (float *)       calloc(nbins, sizeof(float));
  prof->Eay     = (float *)       calloc(nbins, sizeof(float));
  prof->Eaz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ebx     = (float *)       calloc(nbins, sizeof(float));
  prof->Eby     = (float *)       calloc(nbins, sizeof(float));
  prof->Ebz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecx     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecy     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ekin    = (float *)       calloc(nbins, sizeof(float));
  prof->Epot    = (float *)       calloc(nbins, sizeof(float));
  prof->M_gas   = (float *)       calloc(nbins, sizeof(float));
  prof->M_star  = (float *)       calloc(nbins, sizeof(float));
  prof->u_gas   = (float *)       calloc(nbins, sizeof(float));
  prof->Z_gas_sh  = (float *)       calloc(nbins, sizeof(float));
  prof->Z_star_sh = (float *)       calloc(nbins, sizeof(float));
}

/*=================================================================================================
 * free_profile()
 *=================================================================================================*/
void free_profile(halo_profile_t *prof)
{
  free(prof->r);
  free(prof->npart);
  free(prof->M_in_r);
  free(prof->ovdens);
  free(prof->dens);
  free(prof->vcirc);
  free(prof->vesc);
  free(prof->sigv);
  free(prof->Lx);
  free(prof->Ly);
  free(prof->Lz);
  free(prof->b);
  free(prof->c);
  free(prof->Eax);
  free(prof->Eay);
  free(prof->Eaz);
  free(prof->Ebx);
  free(prof->Eby);
  free(prof->Ebz);
  free(prof->Ecx);
  free(prof->Ecy);
  free(prof->Ecz);
  free(prof->Ekin);
  free(prof->Epot);
  free(prof->M_gas);
  free(prof->M_star);
  free(prof->u_gas);
  free(prof->Z_gas_sh);
  free(prof->Z_star_sh);
}

/*=================================================================================================
 * get_TotalNumHalos_from_particles()
 *=================================================================================================*/
uint64_t get_TotalNumHalos_from_particles(int nfiles, char** infile)
{
  uint64_t       FileNumHalos, TotalNumHalos;
  int            ifile;
  int32_t        one;
  int            swap;
  FILE          *fpin;
  
  TotalNumHalos = 0;
  
  for(ifile=0; ifile<nfiles; ifile++) {
    
    // open input file
    fpin = fopen(infile[ifile],"rb");
    if(fpin == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for writing\n",infile[ifile]);
      exit(0);
    }
    
    // figure out swap status
    fread(&one, sizeof(int32_t), 1, fpin);
    if(one == 1)   swap = 0;
    else           swap = 1;
    
    ReadULong(fpin, &FileNumHalos,   swap);
    TotalNumHalos += FileNumHalos;
    
    fclose(fpin);
  }
  
#ifdef VERBOSE
  fprintf(stderr,"o found %ld halos across all files\n",TotalNumHalos);
#endif
  
  return TotalNumHalos;
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned long integer
 */
int ReadULong(FILE *fptr,unsigned long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(unsigned long) == 4)
   {
    if (fread(n,4,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
   }
  else if(sizeof(unsigned long) == 8)
   {
    if (fread(n,8,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
    }
   }
  else
   {
    fprintf(stderr,"ReadULong: something wrong...cannot read long\n");
    exit(0);
   }
  
  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(float) != 4)
   {
    fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap)
   {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
   }
  return(TRUE);
}

