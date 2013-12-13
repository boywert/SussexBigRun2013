
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>


#ifndef LONGIDS
typedef unsigned int MyIDType;
#ifdef HAVE_HDF5
#define HDF5_ID_TYPE H5T_NATIVE_UINT
#endif
#else
typedef unsigned long long MyIDType;
#ifdef HAVE_HDF5
#define HDF5_ID_TYPE H5T_NATIVE_ULLONG
#endif
#endif

extern int SnapFormat;

extern int  SnapSkipFac;
extern int  SnapshotNum;
extern long long TotNumPart;

extern char OutputDir[512];
extern char SnapshotFileBase[512];

extern int  LastSnapShotNr;


#define GENS 2

extern struct halo_catalogue
{
      long long TotNids;
      MyIDType *SubOffset;
      MyIDType *IdList;
      int TotNsubhalos;
      int TotNgroups;
      int *SubLen;
      int *SubParentHalo;
      int *IdToHalo;
      struct descendant_data
      {
           int SnapNum[GENS];
           int HaloIndex[GENS];
#ifdef SKIP_BY_WEIGHT
           float Weight[GENS];
#endif
      }    *Descendant;
      int *CountProgenitors;
}
CatA, CatB, CatC;




extern struct io_header
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions. 
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs 
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.   
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

  char fill[48];		/*!< fills to 256 Bytes */
}
header;				/*!< holds header for snapshot files */


struct twoids
{
  MyIDType id, ord;
};





#endif




