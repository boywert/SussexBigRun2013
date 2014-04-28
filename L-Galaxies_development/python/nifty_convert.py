import pylab
import numpy
import read_lgal
import LGalaxyStruct

def convert_nifty():
    boxsize = 62.5*1000. #kpc/h
    Mpc2kpc = 1000.0
    Gadget2Msun = 1.e10

    folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    

    firstfile = 0
    lastfile = 0
    # f = open(snaplist_file)
    # lines = f.readlines()
    # f.close()
    # i=0
    filter = LGalaxyStruct.properties_used

    filter['Type'] = True
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
    for galaxy in gal:
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
        Type = galaxy["Type"]

        if(Type != 2):
            lasthaloid = haloid
        else:
            last_main_prog = galaxy["FirstProgGal"]
            while (gal[last_main_prog]["Type"] == 2 && last_main_prog > -1)
                last_main_prog = gal[last_main_prog]["FirstProgGal"]
            lasthaloid = gal[last_main_prog]["HaloID"]
        # if(haloid == 60000000000176):
        #     print haloid,Type, Mstar, X, Y, Z, VX, VY, VZ, Mcold, Mhot, Mbh, Z_gas, Z_stars, T_stars

convert_nifty()
