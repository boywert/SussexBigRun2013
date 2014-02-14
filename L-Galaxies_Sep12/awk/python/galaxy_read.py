import LGalaxyStruct
import numpy


folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
file_prefix = "SA_z7.66"
firstfile = 0
lastfile = 0
 
filename = folder+file_prefix+"_"+"%d"%(firstfile)
f = open(filename,"rb")
dummy = numpy.fromfile(f,numpy.int32,1)
nTrees = dummy[0]
dummy = numpy.fromfile(f,numpy.int32,1)
nHalos = dummy[0]
nTreeHalos = numpy.fromfile(f,numpy.int32,nTrees)
Galaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,nHalos)
#print LGalaxyStruct.struct_dtype
f.close()

properties = LGalaxyStruct.properties_used
properties['DiskMass'] = True
out = numpy.dtype([])

for el in LGalaxyStruct.properties_used.keys():
    if(LGalaxyStruct.properties_used[el] is True):
        print el, LGalaxyStruct.properties_used[el]
        numpy.append(out,(el,LGalaxyStruct.struct_dtype[el]))
    
print out
