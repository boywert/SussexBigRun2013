struct LGalaxy {
   int Type;
   int HaloIndex;
   int SnapNum;
   float LookBackTimeToSnap;
   float CentralMvir;
   float Pos[3];
   float Vel[3];
   int Len;
   float Mvir;
   float Rvir;
   float Vvir;
   float Vmax;
   float GasSpin[3];
   float StellarSpin[3];
   float InfallVmax;
   int InfallSnap;
   float HotRadius;
   float OriMergTime;
   float MergTime;
   float DistanceToCentralGal[3];
   float ColdGas;
   float BulgeMass;
   float DiskMass;
   float HotGas;
   float EjectedMass;
   float BlackHoleMass;
   float BlackHoleGas;
   float ICM;
   float MetalsColdGas;
   float MetalsBulgeMass;
   float MetalsDiskMass;
   float MetalsHotGas;
   float MetalsEjectedMass;
   float MetalsICM;
   float Sfr;
   float SfrBulge;
   float XrayLum;
   float BulgeSize;
   float StellarDiskRadius;
   float GasDiskRadius;
   float CosInclination;
   int DisruptOn;
   int MergeOn;
   float CoolingRadius;
   float QuasarAccretionRate;
   float RadioAccretionRate;
   float Mag[15];
   float MagBulge[15];
   float MagDust[15];
   float MassWeightAge;
   int sfh_ibin;
   float sfh_time[19];
   float sfh_dt[19];
   float sfh_DiskMass[19];
   float sfh_BulgeMass[19];
   float sfh_ICM[19];
   float sfh_MetalsDiskMass[19];
   float sfh_MetalsBulgeMass[19];
   float sfh_MetalsICM[19];
};
struct MoMaFGalaxy {
};
