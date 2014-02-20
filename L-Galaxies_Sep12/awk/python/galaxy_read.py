import LGalaxyStruct
import numpy


folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
snaplist_file = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/cubep3m_zlist_out"
firstfile = 0
lastfile = 0





def readsnap(folder,file_prefix,firstfile,lastfile,filter_arr=LGalaxyStruct.properties_used):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    Galaxy = numpy.array([],dtype=LGalaxyStruct.struct_dtype)
    filter_tuple = []
    for prop in filter_arr.keys():
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,LGalaxyStruct.struct_dtype[prop]))
    filter_dtype = numpy.dtype(filter_tuple)
    print filter_dtype
    # output_Galaxy = numpy.array([],dtype=numpy.dtype([("DiskMass",numpy.float32)]))
    # for ifile in range(firstfile,lastfile+1):
    #     #pr
    #     filename = folder+file_prefix+"_"+"%d"%(ifile)
    #     f = open(filename,"rb")
    #     dummy = numpy.fromfile(f,numpy.int32,1)
    #     this_nTrees =  dummy[0]
    #     nTrees += this_nTrees
    #     dummy = numpy.fromfile(f,numpy.int32,1)
    #     this_nHalos = dummy[0]
    #     nHalos += this_nHalos
    #     addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
    #     nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
    #     addedGalaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,this_nHalos)
    #     StellarMass += numpy.sum(addedGalaxy[:]["DiskMass"]) + numpy.sum(addedGalaxy[:]["BulgeMass"])
    #     Sfr += numpy.sum(addedGalaxy[:]["Sfr"])
    #     BlackHoleMass += numpy.sum(addedGalaxy[:]["BlackHoleMass"])
    #     #this_addedGalaxy = numpy.array(addedGalaxy[:]["DiskMass"],dtype=numpy.dtype([("DiskMass",numpy.float32)]))
    #     #output_Galaxy = numpy.append(output_Galaxy,this_addedGalaxy)
    #     Galaxy = numpy.append(Galaxy,addedGalaxy)
    #     #print nTrees,nHalos,len(Galaxy)
    #     f.close()
    # return Galaxy
    # #return (StellarMass ,BlackHoleMass,Sfr)

f = open(snaplist_file)
lines = f.readlines()
f.close()
i=0

filter = LGalaxyStruct.properties_used
filter['Mvir'] = True
filter['Mag'] = True
for this_line in lines:
    i+=1
    if(i != 1):
        #print "z",this_line.strip()
        file_prefix = "SA_z"+this_line.strip()
        readsnap(folder,file_prefix,firstfile,lastfile,filter)
        #print 1./(float(this_line.strip())+1.), StellarMass ,BlackHoleMass,Sfr
#print LGalaxyStruct.struct_dtype
 
