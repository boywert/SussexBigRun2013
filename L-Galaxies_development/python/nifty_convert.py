import pylab
import numpy
import read_lgal
import LGalaxyStruct

def convert_nifty():
    boxsize = 47.0

    folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"


    firstfile = 0
    lastfile = 0
    # f = open(snaplist_file)
    # lines = f.readlines()
    # f.close()
    # i=0
    filter = LGalaxyStruct.properties_used
    filter['Mvir'] = False
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

    print filter
    file_prefix = "SA_"    
    (nGals,gal) = read_lgal.read_lgaltree(folder,file_prefix,firstfile,lastfile,filter)
    for galaxy in gal:
        haloid = galaxy["HaloID"]
        Mstar = galaxy['DiskMass']+galaxy['BulgeMass']
        X = galaxy["Pos"][0]
        Y = galaxy["Pos"][1]
        Z = galaxy["Pos"][2]
        VX = galaxy["Vel"][0]
        VY = galaxy["Vel"][1]
        VZ = galaxy["Vel"][2]
        Mcold = galaxy["ColdGas"]
        Mhot = galaxy["HotGas"]
        Mbh = galaxy["BlackHoleMass"]
        Z_gas = galaxy["MetalsColdGas"] + galaxy["MetalsHotGas"]
        Z_stars = galaxy["MetalsBulgeMass"] + galaxy["MetalsDiskMass"]
        T_stars = galaxy["MassWeightAge"]
        

convert_nifty()
