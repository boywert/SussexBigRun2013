import numpy
import os

lgal_fmt = []
if NO_PROPS_OUTPUTS:
    if GALAXYTREE:
        lgal_fmt.append(("GalID","i8")) #long long GalID;		// ID of galaxy, unique within simulation and SAM run.

if not NO_PROPS_OUTPUTS:
    if GALAXYTREE:
        lgal_fmt.append(("GalID","i8")) #long long GalID; /** ID of galaxy, unique within simulation and SAM run.*/
        lgal_fmt.append(("HaloID","i8")) #long long HaloID; // Unique ID of MPA halo containing this galaxy
    if MBPID:
        lgal_fmt.append(("MostBoundID","i8")) #long long MostBoundID; // Most bound particle at centre of subhalo last associated with this galaxy.  Put here as want all 8-byte blocks together at top of output record.
    if GALAXYTREE:
        lgal_fmt.append(("FirstProgGal","i8")) #long long FirstProgGal;	// Main progenitor of this galaxy. Also the first progenitor in a linked list representation of the merger tree.
        lgal_fmt.append(("NextProgGal","i8")) #long long NextProgGal;	// Next progenitor of this galaxy in linked list representation of merger tree
        lgal_fmt.append(("GalID","i8")) #long long LastProgGal;	// Galaxies with id between this galaxyId and this lastProgenitorId form the merger tree rooted in this galaxy.
        lgal_fmt.append(("FOFCentralGal","i8")) #long long FOFCentralGal;	// TODO has been coded, but that code must be tested.
        lgal_fmt.append(("FileTreeNr","i8")) #long long FileTreeNr;
        lgal_fmt.append(("DescendantGal","i8")) #long long DescendantGal;	// Pointer to the descendant of this galaxy in its merger tree; -1 if there is no descendant
        lgal_fmt.append(("MainLeafId","i8")) #long long MainLeafId;
        lgal_fmt.append(("TreeRootId","i8")) #long long TreeRootId;
        lgal_fmt.append(("GalID","i8")) #long long SubID;
        lgal_fmt.append(("MMSubID","i8")) #long long MMSubID; // fofId, the subhaloid of the subhalo at the center of the fof group
        lgal_fmt.append(("GalID","i4")) #int   PeanoKey; // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
        lgal_fmt.append(("Redshift","f4")) #float Redshift; // redshift of the snapshot where this galaxy resides
    lgal_fmt.append(("Type","i4")) #int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
    if not GALAXYTREE:
        lgal_fmt.append(("HaloIndex","i4")) #int   HaloIndex;
        #// long long SubID;
        #// long long FirstHaloInFOFgroup;
    if HALOPROPERTIES:
        lgal_fmt.append(("HaloM_Mean200","f4")) #float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
        lgal_fmt.append(("HaloM_Crit200","f4"))
        lgal_fmt.append(("HaloM_TopHat","f4"))
        lgal_fmt.append(("HaloPos","3f4")) #float HaloPos[3];
        lgal_fmt.append(("HaloVel","3f4")) #float HaloVel[3];
        lgal_fmt.append(("HaloVelDisp","f4")) #float HaloVelDisp;
        lgal_fmt.append(("HaloVmax","f4")) #float HaloVmax;
        lgal_fmt.append(("HaloSpin","3f4")) #float HaloSpin[3];
    lgal_fmt.append(("SnapNum","i4")) #int   SnapNum; // The snapshot number where this galaxy was identified.
    lgal_fmt.append(("LookBackTimeToSnap","f4")) #float LookBackTimeToSnap; //The time from a given snapshot to z=0, in years
    lgal_fmt.append(("CentralMvir","f4")) #float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
 
   #/* properties of subhalo at the last time this galaxy was a central galaxy */
    lgal_fmt.append(("Pos","3f4")) #float Pos[3]; // 1/h Mpc - Galaxy Positions
    lgal_fmt.append(("Vel","3f4")) #float Vel[3]; // km/s - Galaxy Velocities
    lgal_fmt.append(("Len","i4")) #int   Len;   
    lgal_fmt.append(("Mvir","f4")) #float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
    lgal_fmt.append(("Rvir","f4")) #float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
    lgal_fmt.append(("Vvir","f4")) #float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
    lgal_fmt.append(("Vmax","f4")) #float Vmax; // km/s - Maximum rotational velocity of the subhalo, or the last value for type 2 galaxies.
    lgal_fmt.append(("GasSpin","3f4")) #float GasSpin[3]; // Gas Spin
    lgal_fmt.append(("StellarSpin","3f4")) #float StellarSpin[3]; // Stellar Spin
    lgal_fmt.append(("InfallVmax","f4")) #float InfallVmax; // km/s - Vmax at infall
    lgal_fmt.append(("InfallSnap","i4")) #int InfallSnap; // Snapnum at infall
    lgal_fmt.append(("HotRadius","f4")) #float HotRadius; //Mpc/h - Radius of the hot gas
    
    #/*dynamical friction merger time*/
    lgal_fmt.append(("OriMergTime","f4")) #float OriMergTime;
    lgal_fmt.append(("MergTime","f4")) #float MergTime;
    lgal_fmt.append(("DistanceToCentralGal","3f4")) #float DistanceToCentralGal[3];
    #/* baryonic reservoirs */
    lgal_fmt.append(("ColdGas","f4")) #float ColdGas; // 10^10/h Msun - Mass in cold gas.
    lgal_fmt.append(("BulgeMass","f4")) #float BulgeMass; // 10^10/h Msun - Mass in the bulge
    lgal_fmt.append(("DiskMass","f4")) #float DiskMass;
    lgal_fmt.append(("HotGas","f4")) #float HotGas; // 10^10/h Msun - Mass in hot gas
    lgal_fmt.append(("EjectedMass","f4")) #float EjectedMass; // 10^10/h Msun - Mass in ejected gas
    lgal_fmt.append(("BlackHoleMass","f4")) #float BlackHoleMass; // 10^10/h Msun - Mass in black hole
    lgal_fmt.append(("StellarSpin","3f4")) #float BlackHoleGas; // 10^10/h Msun - Mass in BH accretion disk
    #/* ICL magnitude and mass*/
    lgal_fmt.append(("ICM","f4")) #float ICM;            // mass in intra-cluster stars, for type 0,1
    if METALS:
        struct metals MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
        struct metals MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
        struct metals MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
        struct metals MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
        struct metals MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
        struct metals MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
    if METALS_SELF:
        struct metals MetalsHotGasSelf; // hot gas metals that come from self
    else:
        lgal_fmt.append(("MetalsColdGas","f4")) #float MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
        lgal_fmt.append(("MetalsBulgeMass","f4")) #float MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
        lgal_fmt.append(("MetalsDiskMass","f4")) #float MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
        lgal_fmt.append(("MetalsHotGas","f4")) #float MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
        lgal_fmt.append(("MetalsEjectedMass","f4")) #float MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
        lgal_fmt.append(("MetalsICM","f4")) #float MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
    if METALS_SELF:
        lgal_fmt.append(("MetalsHotGasSelf","f4")) #float MetalsHotGasSelf; // hot gas metals that come from self
    #/* misc */
    lgal_fmt.append(("Sfr","f4")) #float Sfr;
    lgal_fmt.append(("SfrBulge","f4")) #float SfrBulge;
    lgal_fmt.append(("XrayLum","f4")) #float XrayLum;
    lgal_fmt.append(("BulgeSize","f4")) #float BulgeSize;
    lgal_fmt.append(("StellarDiskRadius","f4")) #float StellarDiskRadius;
    lgal_fmt.append(("GasDiskRadius","f4")) #float GasDiskRadius;
    lgal_fmt.append(("CosInclination","f4")) #float CosInclination; // cos(angle) between galaxy spin and the z-axis
    lgal_fmt.append(("StellarSpin","i4")) #int   DisruptOn; // 0: galaxy merged onto merger center; 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center
    if MERGE01:
        lgal_fmt.append(("MergeOn","i4")) #int   MergeOn;   // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....
    lgal_fmt.append(("CoolingRadius","f4")) #float CoolingRadius;  // Q: store this ? (was stored in Delucia20006a)
    lgal_fmt.append(("StellarSpin","f4")) #float QuasarAccretionRate;
    float RadioAccretionRate;


  /* magnitudes in various bands */
if OUTPUT_REST_MAGS:
    float Mag[NMAG]; // rest-frame absolute mags
    float MagBulge[NMAG]; // rest-frame absolute mags for the bulge
    float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
    float MassWeightAge;
    if POST_PROCESS_MAGS:
        float rbandWeightAge;
    if ICL:
        float MagICL[NMAG];          // rest-frame absolute mags of ICL

if OUTPUT_OBS_MAGS:
    float ObsMag[NMAG]; // obs-frame absolute mags
    float ObsMagBulge[NMAG]; // obs-frame absolute mags for the bulge
    float ObsMagDust[NMAG]; // dust-corrected, obs-frame absolute mags
    if ICL:
        float ObsMagICL[NMAG];  // observer-frame absolute mags for intra-cluster light
    if OUTPUT_MOMAF_INPUTS:
        float dObsMag[NMAG];
        float dObsMagBulge[NMAG];
        float dObsMagDust[NMAG];
        if ICL:
            float dObsMagICL[NMAG];       

if STAR_FORMATION_HISTORY:
    int sfh_ibin; //Index of highest bin currently in use
    float sfh_time[SFH_NBIN]; //time to present at the middle of bin in years.
    float sfh_dt[SFH_NBIN]; //time width of bin in years.
    float sfh_DiskMass[SFH_NBIN];
    float sfh_BulgeMass[SFH_NBIN];
    float sfh_ICM[SFH_NBIN];
    if METALS:
        struct metals sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
        struct metals sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
        struct metals sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
    else:
        float sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
        float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
        float sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.


if YIELDS:
    struct elements sfh_ElementsDiskMass[SFH_NBIN];
    struct elements sfh_ElementsBulgeMass[SFH_NBIN];
    struct elements sfh_ElementsICM[SFH_NBIN];
    
    struct elements DiskMass_elements;
    struct elements BulgeMass_elements;
    struct elements ColdGas_elements;
    struct elements HotGas_elements;
    struct elements ICM_elements;
    struct elements EjectedMass_elements;

