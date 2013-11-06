;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LGalaxy_testfloats, LGs, nstart 
; test whether floats are NaN or too small for SQLServer
; assumes the existence of a function testFloat accepting an array of floats
 badranges = []
 bad = 0
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'LookBackTimeToSnap --- ', nstart+sel
     print, 'LookBackTimeToSnap --- ', LGs[sel].LookBackTimeToSnap
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralMvir --- ', nstart+sel
     print, 'CentralMvir --- ', LGs[sel].CentralMvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[0] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[1] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[2] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[0] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[1] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[2] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mvir --- ', nstart+sel
     print, 'Mvir --- ', LGs[sel].Mvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Rvir --- ', nstart+sel
     print, 'Rvir --- ', LGs[sel].Rvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vvir --- ', nstart+sel
     print, 'Vvir --- ', LGs[sel].Vvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vmax --- ', nstart+sel
     print, 'Vmax --- ', LGs[sel].Vmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[0] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[1] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[2] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[0] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[1] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[2] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmax --- ', nstart+sel
     print, 'InfallVmax --- ', LGs[sel].InfallVmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotRadius --- ', nstart+sel
     print, 'HotRadius --- ', LGs[sel].HotRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'OriMergTime --- ', nstart+sel
     print, 'OriMergTime --- ', LGs[sel].OriMergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MergTime --- ', nstart+sel
     print, 'MergTime --- ', LGs[sel].MergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[0] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[1] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[2] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGas --- ', nstart+sel
     print, 'ColdGas --- ', LGs[sel].ColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeMass --- ', nstart+sel
     print, 'BulgeMass --- ', LGs[sel].BulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskMass --- ', nstart+sel
     print, 'DiskMass --- ', LGs[sel].DiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotGas --- ', nstart+sel
     print, 'HotGas --- ', LGs[sel].HotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'EjectedMass --- ', nstart+sel
     print, 'EjectedMass --- ', LGs[sel].EjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BlackHoleMass --- ', nstart+sel
     print, 'BlackHoleMass --- ', LGs[sel].BlackHoleMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BlackHoleGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BlackHoleGas --- ', nstart+sel
     print, 'BlackHoleGas --- ', LGs[sel].BlackHoleGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ICM --- ', nstart+sel
     print, 'ICM --- ', LGs[sel].ICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsColdGas --- ', nstart+sel
     print, 'MetalsColdGas --- ', LGs[sel].MetalsColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsBulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsBulgeMass --- ', nstart+sel
     print, 'MetalsBulgeMass --- ', LGs[sel].MetalsBulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsDiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsDiskMass --- ', nstart+sel
     print, 'MetalsDiskMass --- ', LGs[sel].MetalsDiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsHotGas --- ', nstart+sel
     print, 'MetalsHotGas --- ', LGs[sel].MetalsHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsEjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsEjectedMass --- ', nstart+sel
     print, 'MetalsEjectedMass --- ', LGs[sel].MetalsEjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsICM --- ', nstart+sel
     print, 'MetalsICM --- ', LGs[sel].MetalsICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BurstMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BurstMass --- ', nstart+sel
     print, 'BurstMass --- ', LGs[sel].BurstMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Sfr --- ', nstart+sel
     print, 'Sfr --- ', LGs[sel].Sfr
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'SfrBulge --- ', nstart+sel
     print, 'SfrBulge --- ', LGs[sel].SfrBulge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'XrayLum --- ', nstart+sel
     print, 'XrayLum --- ', LGs[sel].XrayLum
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeSize --- ', nstart+sel
     print, 'BulgeSize --- ', LGs[sel].BulgeSize
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarDiskRadius --- ', nstart+sel
     print, 'StellarDiskRadius --- ', LGs[sel].StellarDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasDiskRadius --- ', nstart+sel
     print, 'GasDiskRadius --- ', LGs[sel].GasDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CosInclination --- ', nstart+sel
     print, 'CosInclination --- ', LGs[sel].CosInclination
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRadius --- ', nstart+sel
     print, 'CoolingRadius --- ', LGs[sel].CoolingRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'QuasarAccretionRate --- ', nstart+sel
     print, 'QuasarAccretionRate --- ', LGs[sel].QuasarAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'RadioAccretionRate --- ', nstart+sel
     print, 'RadioAccretionRate --- ', LGs[sel].RadioAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[0] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[1] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[2] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[3] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[4] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[0] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[1] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[2] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[3] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[4] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[0] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[1] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[2] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[3] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[4] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassWeightAge --- ', nstart+sel
     print, 'MassWeightAge --- ', LGs[sel].MassWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_time(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_time[0] --- ', nstart+sel
     print, 'sfh_time --- ', LGs[sel].sfh_time(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_time(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_time[1] --- ', nstart+sel
     print, 'sfh_time --- ', LGs[sel].sfh_time(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_time(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_time[2] --- ', nstart+sel
     print, 'sfh_time --- ', LGs[sel].sfh_time(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_dt(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_dt[0] --- ', nstart+sel
     print, 'sfh_dt --- ', LGs[sel].sfh_dt(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_dt(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_dt[1] --- ', nstart+sel
     print, 'sfh_dt --- ', LGs[sel].sfh_dt(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_dt(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_dt[2] --- ', nstart+sel
     print, 'sfh_dt --- ', LGs[sel].sfh_dt(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[0] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[1] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[2] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[0] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[1] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[2] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[0] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[1] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[2] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[0] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[1] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[2] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[0] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[1] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[2] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[0] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[1] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[2] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BurstMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BurstMass[0] --- ', nstart+sel
     print, 'sfh_BurstMass --- ', LGs[sel].sfh_BurstMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BurstMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BurstMass[1] --- ', nstart+sel
     print, 'sfh_BurstMass --- ', LGs[sel].sfh_BurstMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BurstMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BurstMass[2] --- ', nstart+sel
     print, 'sfh_BurstMass --- ', LGs[sel].sfh_BurstMass(2)
     badranges=[badranges,sel]
 endif
if(bad) then begin 
     print, 'badranges found: ',badranges
endif
return, badranges
end
