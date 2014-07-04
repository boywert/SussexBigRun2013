import LGalaxyStruct
import numpy
import os

struct_lgalinput = numpy.dtype([
('Descendant',numpy.int32,1),
('FirstProgenitor',numpy.int32,1),
('NextProgenitor',numpy.int32,1),
('FirstHaloInFOFgroup',numpy.int32,1),
('NextHaloInFOFgroup',numpy.int32,1),
('Len',numpy.int32,1),
('M_Mean200',numpy.float32,1),
('M_Crit200',numpy.float32,1),
('M_TopHat',numpy.float32,1),
('Pos',numpy.float32,3),
('Vel',numpy.float32,3),
('VelDisp',numpy.float32,1),
('Vmax',numpy.float32,1),
('Spin',numpy.float32,3),
('MostBoundID',numpy.int64,1),
('SnapNum',numpy.int32,1),
('FileNr',numpy.int32,1),
('SubhaloIndex',numpy.int32,1),
('SubHalfMass',numpy.int32,1)
])

def read_lgalinput(folder,firstfile,lastfile,lastsnap):
    nHalos = 0
    nTrees = 0
    ngalstree = numpy.array([],dtype=numpy.int32)
    output_trees = numpy.array([],dtype=struct_lgalinput)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+"/trees_%03d.%d"%(lastsnap,ifile)
        f = open(filename,"rb")
        this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
        this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
        this_ngalstree = numpy.fromfile(f,numpy.int32,this_nTrees)
        this_trees = numpy.fromfile(f,struct_lgalinput,this_nHalos)
        nHalos += this_nHalos
        nTrees += this_nTrees
        ngalstree = numpy.append(ngalstree,this_ngalstree)
        output_trees = numpy.append(output_trees,this_trees)

    return (nHalos,nTrees,ngalstree,output_trees)



# This function return (nTrees,nHalos,nTreeHalos,Galaxy)
# The input are (folder,file_prefix,firstfile,lastfile [,filter_arr])
def read_lgaltree(folder,file_prefix,firstfile,lastfile,filter_arr=LGalaxyStruct.properties_used):
    nHalos = 0
    filter_tuple = []
    for prop in LGalaxyStruct.struct_dtype.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,LGalaxyStruct.struct_dtype[prop]))
    filter_dtype = numpy.dtype(filter_tuple)

    output_Galaxy = numpy.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"galtree_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        one = dummy[0]
        dummy = numpy.fromfile(f,numpy.int32,1)
        structsize = dummy[0]
        if(structsize != LGalaxyStruct.struct_dtype.itemsize):
            print "size mismatch:",structsize,LGalaxyStruct.struct_dtype.itemsize
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        f.seek(structsize, os.SEEK_SET) 
        print "File ", ifile," nGals = ",this_nHalos
        this_addedGalaxy = numpy.fromfile(f,LGalaxyStruct.struct_dtype,this_nHalos)
        addedGalaxy = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in LGalaxyStruct.struct_dtype.names:
            if(filter_arr[prop] is True):
                addedGalaxy[prop] = this_addedGalaxy[prop]
        output_Galaxy = numpy.append(output_Galaxy,addedGalaxy)
       
      
        f.close()

    print output_Galaxy["Type"]
    return (nHalos,output_Galaxy)


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
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nTrees =  dummy[0]
        nTrees += this_nTrees
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        print "File ", ifile," nGals = ",this_nHalos
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
      
 
