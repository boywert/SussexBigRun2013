import numpy
import cosmocalc
import read_lgal
import os
import sys

folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
snap = "/mnt/lustre/scratch/cs390/47Mpc/snap.txt"

if(len(sys.argv) < 2):
    os.exit()

this_thread = sys.argv[1]

firstfile = this_thread
lastfile = this_thread
lastsnap = 75
h = 0.7
mpart = 0.000150777086098*1.e10/h  #Msun 

a = numpy.loadtxt(snap)
redshift_list = 1./a - 1.
age_list = []
f = []
i = 0
for z in redshift_list:
    z_str="%.3f"%(z)
    output_folder = "/mnt/lustre/scratch/cs390/47Mpc/dmdt_new/"+z_str
    os.system("mkdir -p "+output_folder)
    f.append(open(output_folder+'/'+str(this_thread),"w+"))
    age_list.append(cosmocalc.cosmocalc(z, H0=h*100., WM=0.27, WV=0.73)['zage_Gyr']*1000.) #Myr
    print cosmocalc.cosmocalc(z, H0=h*100., WM=0.27, WV=0.73)
    i += 1



struct_dmdt = numpy.dtype([
    ('Mass',numpy.float32,1),
    ('dmdt',numpy.float32,1)
])


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
                age_cur = age_list[this_treedata[igal]['SnapNum']]
                # nexthalo = this_treedata[igal]['NextHaloInFOFgroup'] 
                # while nexthalo > -1:
                #     print nexthalo,ngalstree[itree],'/',len(this_treedata)
                #    TotalLen += this_treedata[nexthalo]['Len']
                #    nexthalo = this_treedata[nexthalo]['NextHaloInFOFgroup'] 
                
                # print float(this_treedata[igal]['Len'])/TotalLen
                if(this_treedata[igal]['FirstProgenitor'] > -1):
                    mass_prog = this_treedata[this_treedata[igal]['FirstProgenitor']]['Len']*mpart
                    age_prog = age_list[this_treedata[this_treedata[igal]['FirstProgenitor']]['SnapNum']]
                    dm = mass_cur - mass_prog  #Msun
                    dt = (age_cur-age_prog)  #Myr
                else:
                    dm = mass_cur
                    dt = age_cur- age_list[this_treedata[igal]['SnapNum']-1]
                
                f[snap].write("%g  %g  %g  %g  %g\n"%(this_treedata[igal]['Pos'][0],this_treedata[igal]['Pos'][1],this_treedata[igal]['Pos'][2],mass_cur,dm/dt))


#print treedata
