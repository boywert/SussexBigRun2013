import pylab
import read_lgal
import LGalaxyStruct

def uv_l_z8():
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

    file_prefix = "SA_z8.06"    
    (nTrees,nHalos,nTreeHalos,gal) = read_lgal.readsnap_lgal(folder,file_prefix,firstfile,lastfile,filter)
    lgal_hist = pylab.histogram(gal['MagDust'][:,5],bins=50,range=(-21.5,-10.))
    
    lgal_hist_y = lgal_hist[0]
    lgal_hist_x = []
    for i in range(len(lgal_hist_y)):
        lgal_hist_x.append((lgal_hist[1][i]+lgal_hist[1][i+1])/2.)
    print lgal_hist_y
    print lgal_hist_x