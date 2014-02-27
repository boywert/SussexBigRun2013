import pylab
import numpy
import read_lgal
import LGalaxyStruct

def uv_l_z8():
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    model2_folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
    nore_folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs_nore/"
    snaplist_file = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/cubep3m_zlist_out"
    observe_folder="/mnt/lustre/scratch/cs390/codes/cubepm_131212_6_1728_47Mpc_ext2/observed_UVL/"
    firstfile = 0
    lastfile = 215
    f = open(snaplist_file)
    lines = f.readlines()
    f.close()
    i=0
    filter = LGalaxyStruct.properties_used
    filter['Mvir'] = False
    filter['Mag'] = True
    filter['MagDust'] = True

    file_prefix = "SA_z8.06"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(model2_folder,file_prefix,firstfile,lastfile,filter)
    lgal_hist_metal = pylab.histogram(gal['MagDust'][:,5],bins=nbins,range=(min_mag,max_mag))
    
    lgal_hist_metal_y = lgal_hist_metal[0]
    lgal_hist_metal_x = []
    for i in range(len(lgal_hist_metal_y)):
        lgal_hist_metal_x.append((lgal_hist_metal[1][i]+lgal_hist_metal[1][i+1])/2.)
    
    lgal_hist_metal_x = lgal_hist_metal_x - 5.*numpy.log10(hubble_h)
    lgal_hist_metal_y = lgal_hist_metal_y / boxsize**3. / (lgal_hist_metal_x[1]-lgal_hist_metal_x[0])
    lgal_hist_total =  pylab.histogram(gal['Mag'][:,5],bins=nbins,range=(min_mag,max_mag))
    lgal_hist_total_y = lgal_hist_total[0]/ boxsize**3. / (lgal_hist_metal_x[1]-lgal_hist_metal_x[0])

    (nTrees,nGals,nTreeGals,gal2) = read_lgal.readsnap_lgal(nore_folder,file_prefix,firstfile,lastfile,filter)
    lgal_hist_metal_nr = pylab.histogram(gal2['MagDust'][:,5],bins=nbins,range=(min_mag,max_mag))
    
    lgal_hist_metal_nr_y = lgal_hist_metal_nr[0]
    lgal_hist_metal_nr_x = []
    for i in range(len(lgal_hist_metal_nr_y)):
        lgal_hist_metal_nr_x.append((lgal_hist_metal_nr[1][i]+lgal_hist_metal_nr[1][i+1])/2.)
    
    lgal_hist_metal_nr_x = lgal_hist_metal_nr_x - 5.*numpy.log10(hubble_h)
    lgal_hist_metal_nr_y = lgal_hist_metal_nr_y / boxsize**3. / (lgal_hist_metal_nr_x[1]-lgal_hist_metal_nr_x[0])
    lgal_hist_total_nr =  pylab.histogram(gal2['Mag'][:,5],bins=nbins,range=(min_mag,max_mag))
    lgal_hist_total_nr_y = lgal_hist_total_nr[0]/ boxsize**3. / (lgal_hist_metal_nr_x[1]-lgal_hist_metal_nr_x[0])


    pylab.rc('text', usetex=True)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.plot(lgal_hist_metal_x,lgal_hist_metal_y,'r-')
    ax.plot(lgal_hist_metal_x,lgal_hist_total_y,'b-')
    ax.plot(lgal_hist_metal_nr_x,lgal_hist_metal_nr_y,'r--')
    ax.plot(lgal_hist_metal_nr_x,lgal_hist_total_nr_y,'b--')
    bouwens2011_file = observe_folder+"bouwens2011_z8.txt"
    bouwens2011 = numpy.loadtxt(bouwens2011_file)
    bouwens2011_x = bouwens2011[:,0]-5.*numpy.log10(hubble_h)
    bouwens2011_y = (10.**bouwens2011[:,1])/hubble_h**3.

    bouwens2011_errorup = (10.**(bouwens2011[:,1] + bouwens2011[:,5]) - 10.**bouwens2011[:,1])/hubble_h**3.
    bouwens2011_errordown = (10.**bouwens2011[:,1] - 10.**(bouwens2011[:,1] + bouwens2011[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2011_x,bouwens2011_y,yerr=bouwens2011_y/10., fmt='o')
    ax.set_yscale("log")
    print nGals, len(gal)
    pylab.show()


def galaxy_stellar_massftn():
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    model2_folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs/"
    nore_folder = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/outputs_nore/"
    snaplist_file = "/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131212_6_1728_47Mpc_ext2/mergertrees/cubep3m_zlist_out"
    observe_folder="/mnt/lustre/scratch/cs390/codes/cubepm_131212_6_1728_47Mpc_ext2/observed_UVL/"
    firstfile = 0
    lastfile = 215

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    file_prefix = "SA_z8.06"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(model2_folder,file_prefix,firstfile,lastfile,filter)
    stellarmass= pylab.histogram(numpy.log10(gal['DiskMass'][:]+gal['BulgeMass'][:]))
    print gal['DiskMass']
    print gal['BulgeMass']
    # massftn_y = stellarmass[0]
    # massftn_x = []
    # for i in range(len(stellarmass[0])):
    #     massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)

    # pylab.rc('text', usetex=True)
    # fig = pylab.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(massftn_x,massftn_y,'r-')
    # file_prefix = "SA_z8.06"    
    # (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(nore_folder,file_prefix,firstfile,lastfile,filter)
    # stellarmass= pylab.histogram(numpy.log10(gal['DiskMass'][:]+gal['BulgeMass'][:]))
    # massftn_y = stellarmass[0]
    # massftn_x = []
    # for i in range(len(stellarmass[0])):
    #     massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    # ax.plot(massftn_x,massftn_y,'b-')
    # ax.set_yscale("log")
    # pylab.show()
