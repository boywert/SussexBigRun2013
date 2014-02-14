import LGalaxyStruct
import numpy


folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
firstfile = 0
lastfile = 215


file_prefix = "SA_z7.66"
StellarMass = 0.0
Sfr = 0.0
BlackHoleMass = 0.0

def readsnap(folder,file_prefix,firstfile,lastfile):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    Galaxy = numpy.array([],dtype=LGalaxyStruct.struct_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        print "Tree",dummy
        nTrees += dummy[0]
        dummy = numpy.fromfile(f,numpy.int32,1)
        nHalos += dummy[0]
        addednTreeHalos = numpy.fromfile(f,numpy.int32,nTrees)
        nTreeHalos = numpy.concatenate((nTreeHalos,addednTreeHalos))
        addedGalaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,nHalos)
        Galaxy = numpy.concatenate((Galaxy,addedGalaxy))
        f.close()
    return (nTrees,nHalos,nTreeHalos,Galaxy)


(nTrees,nHalos,nTreeHalos,Galaxy) = readsnap(folder,file_prefix,firstfile,lastfile)
StellarMass += numpy.sum(Galaxy[:]["BulgeMass"]) + numpy.sum(Galaxy[:]["DiskMass"])
BlackHoleMass += numpy.sum(Galaxy[:]["BlackHoleMass"])
Sfr += numpy.sum(Galaxy[:]["Sfr"])
print StellarMass,BlackHoleMass,Sfr
print LGalaxyStruct.struct_dtype
 
