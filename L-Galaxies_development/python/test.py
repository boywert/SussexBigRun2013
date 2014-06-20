import numpy
import read_lgal
import LGalaxyStruct

gadget2msun=10.e10
boxsize = 47.0
hubble_h = 0.7
filter = LGalaxyStruct.properties_used
filter['DiskMass'] = True


(nTrees1,nGals1,nTreeGals1,gal1) = read_lgal.readsnap_lgal("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","SA_z6.00",10,10,filter)
(nTrees2,nGals2,nTreeGals2,gal2) = read_lgal.readsnap_lgal("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/","SA_z6.00",10,10,filter)
    
for i in range(nGals1):
    print gal1["DiskMass"][i],gal2["DiskMass"][i]
