CREATE TABLE GALAXIES (
 type INTEGER NOT NULL 
,  haloIndex INTEGER NOT NULL 
,  snapNum INTEGER NOT NULL 
,  lookBackTimeToSnap REAL NOT NULL 
,  centralMvir REAL NOT NULL 
,  centralRvir REAL NOT NULL 
,  pos_1 REAL NOT NULL 
,  pos_2 REAL NOT NULL 
,  pos_3 REAL NOT NULL 
,  vel_1 REAL NOT NULL 
,  vel_2 REAL NOT NULL 
,  vel_3 REAL NOT NULL 
,  len INTEGER NOT NULL 
,  mvir REAL NOT NULL 
,  rvir REAL NOT NULL 
,  vvir REAL NOT NULL 
,  vmax REAL NOT NULL 
,  gasSpin_1 REAL NOT NULL 
,  gasSpin_2 REAL NOT NULL 
,  gasSpin_3 REAL NOT NULL 
,  stellarSpin_1 REAL NOT NULL 
,  stellarSpin_2 REAL NOT NULL 
,  stellarSpin_3 REAL NOT NULL 
,  infallVmax REAL NOT NULL 
,  infallSnap INTEGER NOT NULL 
,  infallHotGas REAL NOT NULL 
,  hotRadius REAL NOT NULL 
,  oriMergTime REAL NOT NULL 
,  mergTime REAL NOT NULL 
,  distanceToCentralGal_1 REAL NOT NULL 
,  distanceToCentralGal_2 REAL NOT NULL 
,  distanceToCentralGal_3 REAL NOT NULL 
,  coldGas REAL NOT NULL 
,  bulgeMass REAL NOT NULL 
,  diskMass REAL NOT NULL 
,  hotGas REAL NOT NULL 
,  ejectedMass REAL NOT NULL 
,  blackHoleMass REAL NOT NULL 
,  blackHoleGas REAL NOT NULL 
,  iCM REAL NOT NULL 
,  metalsColdGas REAL NOT NULL 
,  metalsBulgeMass REAL NOT NULL 
,  metalsDiskMass REAL NOT NULL 
,  metalsHotGas REAL NOT NULL 
,  metalsEjectedMass REAL NOT NULL 
,  metalsICM REAL NOT NULL 
,  primordialAccretionRate REAL NOT NULL 
,  coolingRate REAL NOT NULL 
,  coolingRate_beforeAGN REAL NOT NULL 
,  sfr REAL NOT NULL 
,  sfrBulge REAL NOT NULL 
,  xrayLum REAL NOT NULL 
,  bulgeSize REAL NOT NULL 
,  stellarDiskRadius REAL NOT NULL 
,  gasDiskRadius REAL NOT NULL 
,  cosInclination REAL NOT NULL 
,  disruptOn INTEGER NOT NULL 
,  mergeOn INTEGER NOT NULL 
,  coolingRadius REAL NOT NULL 
,  quasarAccretionRate REAL NOT NULL 
,  radioAccretionRate REAL NOT NULL 
,  mag_1 REAL NOT NULL 
,  mag_2 REAL NOT NULL 
,  mag_3 REAL NOT NULL 
,  mag_4 REAL NOT NULL 
,  mag_5 REAL NOT NULL 
,  mag_6 REAL NOT NULL 
,  mag_7 REAL NOT NULL 
,  mag_8 REAL NOT NULL 
,  mag_9 REAL NOT NULL 
,  mag_10 REAL NOT NULL 
,  mag_11 REAL NOT NULL 
,  mag_12 REAL NOT NULL 
,  mag_13 REAL NOT NULL 
,  mag_14 REAL NOT NULL 
,  mag_15 REAL NOT NULL 
,  magBulge_1 REAL NOT NULL 
,  magBulge_2 REAL NOT NULL 
,  magBulge_3 REAL NOT NULL 
,  magBulge_4 REAL NOT NULL 
,  magBulge_5 REAL NOT NULL 
,  magBulge_6 REAL NOT NULL 
,  magBulge_7 REAL NOT NULL 
,  magBulge_8 REAL NOT NULL 
,  magBulge_9 REAL NOT NULL 
,  magBulge_10 REAL NOT NULL 
,  magBulge_11 REAL NOT NULL 
,  magBulge_12 REAL NOT NULL 
,  magBulge_13 REAL NOT NULL 
,  magBulge_14 REAL NOT NULL 
,  magBulge_15 REAL NOT NULL 
,  magDust_1 REAL NOT NULL 
,  magDust_2 REAL NOT NULL 
,  magDust_3 REAL NOT NULL 
,  magDust_4 REAL NOT NULL 
,  magDust_5 REAL NOT NULL 
,  magDust_6 REAL NOT NULL 
,  magDust_7 REAL NOT NULL 
,  magDust_8 REAL NOT NULL 
,  magDust_9 REAL NOT NULL 
,  magDust_10 REAL NOT NULL 
,  magDust_11 REAL NOT NULL 
,  magDust_12 REAL NOT NULL 
,  magDust_13 REAL NOT NULL 
,  magDust_14 REAL NOT NULL 
,  magDust_15 REAL NOT NULL 
,  massWeightAge REAL NOT NULL 
 -- size = 428
)
CREATE TABLE SFH_Times (
 -- size = 0
)
