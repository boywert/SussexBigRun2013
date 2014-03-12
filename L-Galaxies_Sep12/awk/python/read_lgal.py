import LGalaxyStruct
import numpy

# This function return (nTrees,nHalos,nTreeHalos,Galaxy)
# The input are (folder,file_prefix,firstfile,lastfile [,filter_arr])
def readsnap_lgal(folder,file_prefix,firstfile,lastfile,filter_arr=LGalaxyStruct.properties_used):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    filter_tuple = []
    for prop in LGalaxyStruct.struct_dtype.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,LGalaxyStruct.struct_dtype[prop]))
    filter_dtype = numpy.dtype(filter_tuple)
    output_Galaxy = numpy.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nTrees =  dummy[0]
        nTrees += this_nTrees
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        print "File" ifile " nGals = " this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
        this_addedGalaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,this_nHalos)
        addedGalaxy = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in LGalaxyStruct.struct_dtype.names:
            if(filter_arr[prop] is True):
                addedGalaxy[prop] = this_addedGalaxy[prop]
        output_Galaxy = numpy.append(output_Galaxy,addedGalaxy)
       
      
        f.close()
    return (nTrees,nHalos,nTreeHalos,output_Galaxy)

def example():
    folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
    snaplist_file = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/cubep3m_zlist_out"
    firstfile = 0
    lastfile = 215
    f = open(snaplist_file)
    lines = f.readlines()
    f.close()
    i=0
    filter = LGalaxyStruct.properties_used
    filter['Mvir'] = False
    filter['Mag'] = False
    filter['MagDust'] = True
    for this_line in lines:
        i+=1
        if(i == 57):
            file_prefix = "SA_z"+this_line.strip()
            (nTrees,nHalos,nTreeHalos,gal) = readsnap_lgal(folder,file_prefix,firstfile,lastfile,filter)
            print gal['MagDust'][:,5]
            #help(gal)
      
 
