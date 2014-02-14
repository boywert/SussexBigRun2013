import LGalaxyStruct
import numpy


folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
snaplist_file = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/cubep3m_zlist_out"
firstfile = 0
lastfile = 215





def readsnap(folder,file_prefix,firstfile,lastfile):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    Galaxy = numpy.array([],dtype=LGalaxyStruct.struct_dtype)
    output_Galaxy = numpy.array([],dtype=numpy.dtype([("DiskMass",numpy.float32)]))
    StellarMass = 0.0
    Sfr = 0.0
    BlackHoleMass = 0.0
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
        StellarMass += numpy.sum(addedGalaxy[:]["DiskMass"]) + numpy.sum(addedGalaxy[:]["BulgeMass"])
        Sfr += numpy.sum(addedGalaxy[:]["Sfr"])
        BlackHoleMass += numpy.sum(addedGalaxy[:]["BlackHoleMass"])
        #this_addedGalaxy = numpy.array(addedGalaxy[:]["DiskMass"],dtype=numpy.dtype([("DiskMass",numpy.float32)]))
        #output_Galaxy = numpy.append(output_Galaxy,this_addedGalaxy)
        #Galaxy = numpy.append(Galaxy,addedGalaxy)
        #print nTrees,nHalos,len(Galaxy)
        f.close()
    return (StellarMass ,BlackHoleMass,Sfr)

f = open(snaplist_file)
lines = f.readlines()
f.close()
i=0
for this_line in lines:
    i+=1
    if(i != 1):
        print "z",this_line.strip()
        file_prefix = "SA_z"+this_line.strip()
        (StellarMass ,BlackHoleMass,Sfr) = readsnap(folder,file_prefix,firstfile,lastfile)
        print StellarMass,BlackHoleMass,Sfr
#print LGalaxyStruct.struct_dtype
 
