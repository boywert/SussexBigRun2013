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
    filter['Mag'] = True
    filter['MagDust'] = True
    print filter
    file_prefix = "SA_"    
    (nGals,gal) = read_lgal.read_lgaltree(folder,file_prefix,firstfile,lastfile,filter)
    for galaxy in gal:
        print galaxy['DiskMass']
    


convert_nifty()
