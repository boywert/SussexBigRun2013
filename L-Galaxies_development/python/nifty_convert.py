import pylab
import numpy
import read_lgal
import LGalaxyStruct
import os

def read_galaxies():


    folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    

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
    filter['ColdGas'] = True
    filter['HotGas'] = True
    filter['HaloID'] = True
    filter['Pos'] = True
    filter['Vel'] = True
    filter['ColdGas'] = True
    filter['HotGas'] = True
    filter['BlackHoleMass'] = True
    filter['MetalsColdGas'] = True
    filter['MetalsHotGas'] = True
    filter['MetalsBulgeMass'] = True
    filter['MetalsDiskMass'] = True
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


AHFdir = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/AHF/"
AHFprefix = "62.5_dm"
SUSSINGtree = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/MergerTree+AHF.txt"
SNAPfile = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/snapidzred.txt"
FileOut = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/trees_061.0"
FileOut2 = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/tree_dbids_061.0"
timesnap = numpy.loadtxt(SNAPfile)

boxsize = 62.5*1000. #kpc/h
Mpc2kpc = 1000.0
Gadget2Msun = 1.e10

for time in timesnap:
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
        for halo in data:
            hid = long(halo[0])
            if hid in halomap:
                print hid, len(halomap[hid])
                for galid in halomap[hid]:
                    galaxy = gal[galid]
                    haloid = galaxy["HaloID"]
                    Mstar = galaxy['DiskMass']+galaxy['BulgeMass']
                    X = galaxy["Pos"][0]*Mpc2kpc
                    Y = galaxy["Pos"][1]*Mpc2kpc
                    Z = galaxy["Pos"][2]*Mpc2kpc
                    VX = galaxy["Vel"][0]
                    VY = galaxy["Vel"][1]
                    VZ = galaxy["Vel"][2]
                    Mcold = galaxy["ColdGas"]*Gadget2Msun
                    Mhot = galaxy["HotGas"]*Gadget2Msun
                    Mbh = galaxy["BlackHoleMass"]*Gadget2Msun
                    Z_gas = (galaxy["MetalsColdGas"] + galaxy["MetalsHotGas"])*Gadget2Msun
                    Z_stars = (galaxy["MetalsBulgeMass"] + galaxy["MetalsDiskMass"])*Gadget2Msun
                    T_stars = galaxy["MassWeightAge"]
                    lasthaloid = lasthalomap[galid]
                    if(galaxy["Type"] == 2):
                        orphan_flag = 1
                    else:
                        orphan_flag = 0
                    print lasthaloid,orphan_flag,X,Y,Z,VX,VY,VZ,Mcold,Mhot,Mbh,Z_gas,Z_stars,T_stars
            else:
                print hid, 0

