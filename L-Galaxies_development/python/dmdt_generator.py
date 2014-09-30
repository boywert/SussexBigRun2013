import numpy
import cosmocalc
import read_lgal

folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
snap = "/mnt/lustre/scratch/cs390/47Mpc/snap.txt"

a = numpy.loadtxt(snap)
redshift_list = 1./a - 1.

firstfile = 127
lastfile = 127
lastsnap = 75
h = 0.7
mpart = 0.000150777086098*1.e10/h  #Msun 

(nHalos,nTrees,ngalstree,treedata) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)

lasthalointree = numpy.cumsum(ngalstree)
firsthalointree = lasthalointree-ngalstree 

#for firsthalo in firsthalointree:
#    firsthalo = int(firsthalo)
#    print treedata[firsthalo]


for itree in range(nTrees):
    this_treedata = numpy.array(treedata[firsthalointree[itree]:lasthalointree[itree]],dtype=read_lgal.struct_lgalinput)
    for snap in range(lastsnap+1):
        for igal in range(ngalstree[itree]):
            if(this_treedata[igal]['SnapNum'] == snap):
                mass_cur = this_treedata[igal]['Len']*mpart
                z_cur = redshift_list[this_treedata[igal]['SnapNum']]
                age_cur = cosmocalc(z_cur, H0=h*100., WM=0.27, WV=0.73)['zage']
                # nexthalo = this_treedata[igal]['NextHaloInFOFgroup'] 
                # while nexthalo > -1:
                #     print nexthalo,ngalstree[itree],'/',len(this_treedata)
                #    TotalLen += this_treedata[nexthalo]['Len']
                #    nexthalo = this_treedata[nexthalo]['NextHaloInFOFgroup'] 
               
                # print float(this_treedata[igal]['Len'])/TotalLen
                mass_prog = this_treedata[this_treedata[igal]['FirstProgenitor']]['Len']*mpart
                z_prog = redshift_list[this_treedata[redshift_list[this_treedata[igal]['SnapNum']]]['SnapNum']]
                age_prog = cosmocalc(z_prog, H0=h*100., WM=0.27, WV=0.73)['zage']
                dm = mass_cur - mass_prog  #Msun
                dt = age_cur-age_prog  #Gyr

                print dm/dt


#print treedata
