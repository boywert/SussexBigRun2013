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
    output_Galaxy = numpy.array([],dtype=numpy.dtype([("DiskMass",numpy.float32])))
    for ifile in range(firstfile,lastfile+1):
        print "File:",ifile
        filename = folder+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nTrees =  dummy[0]
        nTrees += this_nTrees
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
        addedGalaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,this_nHalos)
        output_addedGalaxy[:]["DiskMass"] = addedGalaxy[:]["DiskMass"]
        output_Galaxy = numpy.append(output_Galaxy,output_addedGalaxy)
        #Galaxy = numpy.append(Galaxy,addedGalaxy)
        #print nTrees,nHalos,len(Galaxy)
        f.close()
    return (nTrees,nHalos,nTreeHalos,output_Galaxy)


(nTrees,nHalos,nTreeHalos,Galaxy) = readsnap(folder,file_prefix,firstfile,lastfile)
StellarMass += numpy.sum(Galaxy[:]["DiskMass"])
#BlackHoleMass += numpy.sum(Galaxy[:]["BlackHoleMass"])
#Sfr += numpy.sum(Galaxy[:]["Sfr"])
print StellarMass #,BlackHoleMass,Sfr
#print LGalaxyStruct.struct_dtype
 
