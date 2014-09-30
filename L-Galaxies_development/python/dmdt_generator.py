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
for itree in range(nTrees):
    this_treedata = numpy.array(treedata[firsthalointree[itree]:firsthalointree[itree]+ngalstree[itree]],dtype=read_lgal.struct_lgalinput)
    for snap in range(lastsnap+1):
        for igal in range(ngalstree[itree]):
            if(this_treedata[igal]['FirstHaloInFOFgroup']== igal) & (this_treedata[igal]['SnapNum'] == snap):
                TotalLen = this_treedata[igal]['Len']
                nexthalo = this_treedata[igal]['NextHaloInFOFgroup'] 
                while nexthalo > -1:
                    TotalLen += this_treedata[nexthalo]['Len']
                    nexthalo = this_treedata[nexthalo]['NextHaloInFOFgroup'] 
                print TotalLen

#print treedata
