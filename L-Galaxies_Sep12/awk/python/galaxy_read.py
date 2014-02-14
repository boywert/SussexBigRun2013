import LGalaxyStruct
import numpy


folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
file_prefix = "SA_z7.66"
firstfile = 0
lastfile = 0
StellarMass = 0.0
for ifile in range(firstfile,lastfile):
    filename = folder+file_prefix+"_"+"%d"%(ifile)
    f = open(filename,"rb")
    dummy = numpy.fromfile(f,numpy.int32,1)
    nTrees = dummy[0]
    dummy = numpy.fromfile(f,numpy.int32,1)
    nHalos = dummy[0]
    nTreeHalos = numpy.fromfile(f,numpy.int32,nTrees)
    Galaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,nHalos)
    f.close()
    print Galaxy

 
