import pylab
import numpy
import sys
sys.path.insert(0, "../")
import read_lgal
import LGalaxyStruct
import os
import config
import time


AHFdir = config.AHFdir
AHFprefix = config.AHFprefix
SUSSINGtree = config.SUSSINGtree
SNAPfile = config.SNAPfile
FileOut = config.FileOut
FileOut2 = config.FileOut2
prefix = config.prefix
output_folder = config.output_folder
lgaltree_output = config.lgaltree_output

def read_galaxies():


    folder = config.lgaltree_output
    

    firstfile = 0
    lastfile = 0
    # f = open(snaplist_file)
    # lines = f.readlines()
    # f.close()
    # i=0
    filter = LGalaxyStruct.properties_used

    filter['Type'] = True
    filter['GalID'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    filter['Mag'] = True
    filter['MagDust'] = True
    filter['ColdGas'] = True
    filter['HotGas'] = True
    filter['HaloID'] = True
    filter['Pos'] = True
    filter['Vel'] = True 
    filter['ColdGas'] = True
    filter['HotGas'] = True
    filter['BlackHoleMass'] = True
    filter['EjectedMass'] = True
    filter['ICM'] = True
    filter['BlackHoleGas'] = True
    filter['MetalsColdGas'] = True
    filter['MetalsHotGas'] = True
    filter['MetalsBulgeMass'] = True
    filter['MetalsDiskMass'] = True
    filter['MetalsEjectedMass'] = True
    filter['MetalsICM'] = True
    filter['Sfr'] = True
    filter['SfrBulge'] = True
    filter['MassWeightAge'] = True
    filter["FirstProgGal"] = True


    file_prefix = "SA_"    
    (nGals,gal) = read_lgal.read_lgaltree(folder,file_prefix,firstfile,lastfile,filter)
    halomap = {}
    lasthalomap = []
    for idgal in range(nGals):
        galaxy = gal[idgal]
        haloid = galaxy["HaloID"]
        Type = galaxy["Type"]
        
        if(Type != 2):
            lasthaloid = haloid
        else:
            last_main_prog = idgal+1
            while (gal[last_main_prog]["Type"] == 2):
                last_main_prog = last_main_prog+1
            lasthaloid = gal[last_main_prog]["HaloID"]
        lasthalomap.append(lasthaloid)

        if haloid in halomap:
            halomap[haloid].append(idgal)
        else:
            halomap[haloid] = []
            halomap[haloid].append(idgal)

    return (nGals,gal,halomap,lasthalomap)


(nGals,gal,halomap,lasthalomap) = read_galaxies()

timesnap = numpy.loadtxt(SNAPfile)

boxsize = 62.5*1000. #kpc/h
Mpc2kpc = 1000.0
Gadget2Msun = 1.e10
L_sun = 3.846e26 #W
L_vega = 40.12*L_sun
pc2m = 3.08567758e16 #m
Zsolar = 0.02

os.system("mkdir -p "+output_folder)
for timeid in range(len(timesnap)):
    time = timesnap[timeid]
    ofilename = output_folder+"/"+prefix+"."+"%04d.txt"%timeid
    fp = open(ofilename,"w+")
    zstring = "%.3f" % (time[2])
        #print zstring[len(zstring)-1]
    filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        #print "checking "+filename
    if os.path.isfile(filename) == False:
        zstring = "%.3f" % (time[2]+0.00001)
            #print zstring[len(zstring)-1]
        filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
    if os.path.isfile(filename) == False:
        zstring = "%.3f" % (time[2]-0.00001)
            #print "checking "+filename
            #print zstring[len(zstring)-1]
        filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
    if os.path.isfile(filename) == False:
        print "checking error"+filename
        exit()
    print "Reading "+filename
    stat = os.stat(filename)
        #print stat.st_size
    if(stat.st_size > 384):
        
        data = numpy.loadtxt(filename)
        shape = data.shape
        if (len(shape) == 1) & (shape[0] > 1):
            data_tmp = []
            data_tmp.append(data)
            data = data_tmp

        for halo in data:
            hid = long(halo[0])
            if hid in halomap:
                for galid in halomap[hid]:
                    galaxy = gal[galid]
                    if(galaxy["Type"] == 2):
                        orphan_flag = 1
                    else:
                        orphan_flag = 0
                    lasthaloid = lasthalomap[galid]
                    X = galaxy["Pos"][0]*Mpc2kpc
                    Y = galaxy["Pos"][1]*Mpc2kpc
                    Z = galaxy["Pos"][2]*Mpc2kpc
                    VX = galaxy["Vel"][0]
                    VY = galaxy["Vel"][1]
                    VZ = galaxy["Vel"][2]
                    Mcold = galaxy["ColdGas"]*Gadget2Msun
                    Mhot = galaxy["HotGas"]*Gadget2Msun
                    Mstar = (galaxy['DiskMass']+galaxy['BulgeMass']) * Gadget2Msun
                    Mbh = galaxy["BlackHoleMass"]*Gadget2Msun
                    Z_gas = (galaxy["MetalsColdGas"] + galaxy["MetalsHotGas"]) / (galaxy["ColdGas"] + galaxy["HotGas"])/Zsolar
		    if(Mstar > 0.0):
                    	Z_stars = (galaxy["MetalsBulgeMass"] + galaxy["MetalsDiskMass"]) / (galaxy['DiskMass']+galaxy['BulgeMass'])/Zsolar
		    else:
			Z_stars = 0.0
                    T_stars = galaxy["MassWeightAge"]
                    SFR = galaxy["Sfr"]
                    SFRbulge = galaxy["SfrBulge"]
                    Mhot_halo = Mhot
                    Mcold_halo = -1
                    Mejected = galaxy["EjectedMass"]*Gadget2Msun
                    M_outflow = -1
                    Mgas_disk = Mcold
                    Mgas_spheroid = -1
                    Mstar_disk = galaxy["DiskMass"]*Gadget2Msun
                    Mstar_spheroid = galaxy["BulgeMass"]*Gadget2Msun
                    M_bh = Mbh
                    M_ICstars = galaxy['ICM']*Gadget2Msun
                    M_total = Mstar_spheroid+Mstar_disk+Mhot+Mcold+Mejected+M_ICstars+M_bh
                    MZhot_halo = galaxy["MetalsHotGas"]*Gadget2Msun
                    MZcold_halo = -1
                    MZejected = galaxy["MetalsEjectedMass"]*Gadget2Msun
                    MZ_outflow = -1
                    MZgas_disk = galaxy["MetalsColdGas"]*Gadget2Msun
                    MZgas_spheroid = -1
                    MZstars_disk  = galaxy["MetalsDiskMass"]*Gadget2Msun
                    MZstars_spheroid = galaxy["MetalsBulgeMass"]*Gadget2Msun
                    MZ_bh = -1
                    MZ_ICstars = galaxy["MetalsICM"]
                    BoverT = galaxy["BulgeMass"]/(galaxy['DiskMass']+galaxy['BulgeMass'])
                    r_half = -1
                    r_half_bulge = -1
                    r_half_disk = -1
                   # Magnitude with dust
                    if(galaxy['MagDust'][0] < 98.):
                        nuv_ext =  galaxy['MagDust'][0]
                    else:
                        nuv_ext = 0.
                    if(galaxy['MagDust'][1] < 98.):
                        B_ext = galaxy['MagDust'][1]  
                    else:
                        B_ext = 0.
                    if(galaxy['MagDust'][2] < 98.):
                        V_ext = galaxy['MagDust'][2]   
                    else:
                        V_ext = 0.
                    if(galaxy['MagDust'][3] < 98.):
                        g_ext = galaxy['MagDust'][3]  
                    else:
                        g_ext = 0.
                    if(galaxy['MagDust'][4] < 98.):
                        r_ext = galaxy['MagDust'][4]  
                    else:
                        r_ext = 0.
                    if(galaxy['MagDust'][5] < 98.):
                        K_ext = galaxy['MagDust'][5]  
                    else:
                        K_ext = 0.


                    # Magnitude without dust
                    if(galaxy['Mag'][0] < 98.):
                        nuv =  galaxy['Mag'][0]
                    else:
                        nuv = 0.
                    if(galaxy['Mag'][1] < 98.):
                        B = galaxy['Mag'][1]  
                    else:
                        B = 0.
                    if(galaxy['Mag'][2] < 98.):
                        V = galaxy['Mag'][2]   
                    else:
                        V = 0.
                    if(galaxy['Mag'][3] < 98.):
                        g = galaxy['Mag'][3]  
                    else:
                        g = 0.
                    if(galaxy['Mag'][4] < 98.):
                        r = galaxy['Mag'][4]  
                    else:
                        r = 0.
                    if(galaxy['Mag'][5] < 98.):
                        K = galaxy['Mag'][5]  
                    else:
                        K = 0.                    
                    print >> fp,hid,lasthaloid,orphan_flag,X,Y,Z,VX,VY,VZ,Mcold,Mhot,Mstar,Mbh,Z_gas,Z_stars,T_stars,SFR,SFRbulge,Mhot_halo,Mcold_halo,Mejected,M_outflow,Mgas_disk,Mgas_spheroid,Mstar_disk,Mstar_spheroid,M_bh,M_ICstars,M_total,MZhot_halo,MZcold_halo,MZejected,MZ_outflow,MZgas_disk,MZgas_spheroid,MZstars_disk,MZstars_spheroid,MZ_bh,MZ_ICstars,BoverT,r_half,r_half_bulge,r_half_disk,nuv_ext,B_ext,V_ext,g_ext,r_ext,K_ext,nuv,B,V,g,r,K

    fp.close()
