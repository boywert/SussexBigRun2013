;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_zerofloats, LGs 
; test whether floats are NaN or too small for SQLServer
; if so, set offending values to 0
; assumes the existence of a function testFloat accepting an array of floats
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     LGs[sel].LookBackTimeToSnap = 0
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralMvir = 0
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(0) = 0
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(1) = 0
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(2) = 0
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(0) = 0
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(1) = 0
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(2) = 0
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Mvir = 0
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Rvir = 0
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Vvir = 0
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     LGs[sel].Vmax = 0
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(0) = 0
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(1) = 0
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(2) = 0
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(0) = 0
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(1) = 0
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(2) = 0
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmax = 0
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].HotRadius = 0
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].OriMergTime = 0
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].MergTime = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(0) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(1) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(2) = 0
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGas = 0
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeMass = 0
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskMass = 0
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].HotGas = 0
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].EjectedMass = 0
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BlackHoleMass = 0
 endif
 sel = testFloat(LGs.BlackHoleGas)
 if(sel(0) gt -1) then begin
     LGs[sel].BlackHoleGas = 0
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     LGs[sel].ICM = 0
 endif
 sel = testFloat(LGs.MetalsColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsColdGas = 0
 endif
 sel = testFloat(LGs.MetalsBulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsBulgeMass = 0
 endif
 sel = testFloat(LGs.MetalsDiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsDiskMass = 0
 endif
 sel = testFloat(LGs.MetalsHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsHotGas = 0
 endif
 sel = testFloat(LGs.MetalsEjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsEjectedMass = 0
 endif
 sel = testFloat(LGs.MetalsICM)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsICM = 0
 endif
 sel = testFloat(LGs.BurstMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BurstMass = 0
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     LGs[sel].Sfr = 0
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     LGs[sel].SfrBulge = 0
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     LGs[sel].XrayLum = 0
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeSize = 0
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarDiskRadius = 0
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].GasDiskRadius = 0
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     LGs[sel].CosInclination = 0
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRadius = 0
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].QuasarAccretionRate = 0
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].RadioAccretionRate = 0
 endif
 sel = testFloat(LGs.ObsMag(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(0) = 0
 endif
 sel = testFloat(LGs.ObsMag(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(1) = 0
 endif
 sel = testFloat(LGs.ObsMag(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(2) = 0
 endif
 sel = testFloat(LGs.ObsMag(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(3) = 0
 endif
 sel = testFloat(LGs.ObsMag(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(4) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(0) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(1) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(2) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(3) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(4) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(0) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(1) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(2) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(3) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(4) = 0
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].MassWeightAge = 0
 endif
 sel = testFloat(LGs.sfh_time(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_time(0) = 0
 endif
 sel = testFloat(LGs.sfh_time(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_time(1) = 0
 endif
 sel = testFloat(LGs.sfh_time(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_time(2) = 0
 endif
 sel = testFloat(LGs.sfh_dt(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_dt(0) = 0
 endif
 sel = testFloat(LGs.sfh_dt(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_dt(1) = 0
 endif
 sel = testFloat(LGs.sfh_dt(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_dt(2) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_BurstMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BurstMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_BurstMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BurstMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_BurstMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BurstMass(2) = 0
 endif
end
