import numpy
import read_lgal

folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
firstfile = 127
lastfile = 127
lastsnap = 75

(nHalos,nTrees,ngalstree,treedata) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)

firsthalointree = numpy.cumsum(ngalstree) - numpy.ones(nTrees)*ngalstree[0]

#for firsthalo in firsthalointree:
#    firsthalo = int(firsthalo)
#    print treedata[firsthalo]
for snap in range(lastsnap+1):
    for itree in range(nTrees):
        this_treedata = numpy.array([],dtype=read_lgal.struct_lgalinput)
        this_treedata = numpy.append(this_treedata,treedata[firsthalointree[itree]:firsthalointree[itree]+ngalstree[itree]])
        for igal in range(ngalstree[itree]):
            if(this_treedata[igal].FirstHaloInFOFgroup == igal) & (this_treedata[igal].SnapNum == snap):
                print this_treedata[igal]

#print treedata
