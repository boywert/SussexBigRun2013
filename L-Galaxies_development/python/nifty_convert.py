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
L_sun = 3.846e26 #W
L_vega = 40.12*L_sun
pc2m = 3.08567758e16 #m
output_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/SAM/"
os.system("mkdir -p "+output_folder)
prefix = "nifty_I"
for timeid in range(len(timesnap)):
    time = timesnap[timeid]
    ofilename = output_folder+"/"+prefix+"."+"%04d.txt"%timeid
    fp = open(ofilename,"w+")
    print >> fp, "M_sun/h  L_Vega"
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
                print >> fp,  hid, len(halomap[hid])
                for galid in halomap[hid]:
                    galaxy = gal[galid]
                    haloid = galaxy["HaloID"]
                    Mstar = (galaxy['DiskMass']+galaxy['BulgeMass']) * Gadget2Msun
                    X = galaxy["Pos"][0]*Mpc2kpc
                    Y = galaxy["Pos"][1]*Mpc2kpc
                    Z = galaxy["Pos"][2]*Mpc2kpc
                    VX = galaxy["Vel"][0]
                    VY = galaxy["Vel"][1]
                    VZ = galaxy["Vel"][2]
                    Mcold = galaxy["ColdGas"]*Gadget2Msun
                    Mhot = galaxy["HotGas"]*Gadget2Msun
                    Mbh = galaxy["BlackHoleMass"]*Gadget2Msun
                    Z_gas = (galaxy["MetalsColdGas"] + galaxy["MetalsHotGas"]) / (galaxy["ColdGas"] + galaxy["HotGas"])
                    Z_stars = (galaxy["MetalsBulgeMass"] + galaxy["MetalsDiskMass"]) / (galaxy['DiskMass']+galaxy['BulgeMass'])
                    T_stars = galaxy["MassWeightAge"]
                    lasthaloid = lasthalomap[galid]
                    # convert SDSS magnitude to luminosity
                    # u,g,r,i,z
                    # M = -2.5 log (L_SDSS/ L_vega))

                    if(galaxy['Mag'][0] < 98.):
                        u_l = 10.**(-1.*(galaxy['Mag'][0] + 0.0)/2.5) 
                    else:
                        u_l = 0.
                    if(galaxy['Mag'][1] < 98.):
                        g_l = 10.**(-1.*(galaxy['Mag'][1] + 0.0)/2.5) 
                    else:
                        g_l = 0.
                    if(galaxy['Mag'][2] < 98.):
                        r_l = 10.**(-1.*(galaxy['Mag'][2] + 0.0)/2.5)  
                    else:
                        r_l = 0.
                    if(galaxy['Mag'][3] < 98.):
                        i_l = 10.**(-1.*(galaxy['Mag'][3] + 0.0)/2.5)  
                    else:
                        i_l = 0.
                    if(galaxy['Mag'][4] < 20.):
                        z_l = 10.**(-1.*(galaxy['Mag'][4] + 0.0)/2.5)  
                    else:
                        z_l = 0.

                    if(galaxy["Type"] == 2):
                        orphan_flag = 1
                    else:
                        orphan_flag = 0
                    
                    print >> fp, lasthaloid,orphan_flag,X,Y,Z,VX,VY,VZ,Mcold,Mhot,Mstar, Mbh, Z_gas,Z_stars,T_stars,u_l,g_l,r_l,i_l,z_l
            else:
                print >>fp, hid, 0
    
    fp.close()
