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
    
    model2_folder = "/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/"
    nore_folder = "/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/"
    snaplist_file = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
    observe_folder="/mnt/lustre/scratch/cs390/codes/47Mpc/observed_UVL/"


    firstfile = 0
    lastfile = 127
    f = open(snaplist_file)
    lines = f.readlines()
    f.close()
    i=0
    filter = LGalaxyStruct.properties_used
    filter['Mvir'] = False
    filter['Mag'] = True
    filter['MagDust'] = True

    file_prefix = "SA_z6.00"    
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
    #ax.plot(lgal_hist_metal_x,lgal_hist_metal_y,'r-',label="MagDust(Okamoto)")
    ax.plot(lgal_hist_metal_x,lgal_hist_total_y,'b-',label="Mag(Okamoto)")
    #ax.plot(lgal_hist_metal_nr_x,lgal_hist_metal_nr_y,'r--',label="MagDust(No Reionisation)")
    ax.plot(lgal_hist_metal_nr_x,lgal_hist_total_nr_y,'b--',label="Mag(No Reionisation)")
    ax = add_obs_uv_z6(observe_folder,ax)
    ax.set_yscale("log")
    
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    print nGals, len(gal)
    ax.set_xlabel(r"Mag$_UV$")
    ax.set_ylabel(r"$\phi$ (Galaxies $h^{3}$ Mpc$^{-3}$ mags$^{-1}$)")
    fig.suptitle("GFUV z~8")
    pylab.savefig('uv_l_z8.pdf',bbox_inches='tight')
    pylab.show()


def add_obs_uv_z8(observe_folder,ax):
    bouwens2011_file = observe_folder+"/bouwens2011_z8.txt"
    bouwens2011 = numpy.loadtxt(bouwens2011_file)
    bouwens2011_x = bouwens2011[:,0]-5.*numpy.log10(hubble_h)
    bouwens2011_y = (10.**bouwens2011[:,1])/hubble_h**3.

    bouwens2011_errorup = (10.**(bouwens2011[:,1] + bouwens2011[:,5]) - 10.**bouwens2011[:,1])/hubble_h**3.
    bouwens2011_errordown = (10.**bouwens2011[:,1] - 10.**(bouwens2011[:,1] + bouwens2011[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2011_x,bouwens2011_y,yerr=bouwens2011_y/10., fmt='o',label="Bouwens et al. (2011)") 
    return ax

def add_obs_uv_z6(observe_folder,ax):
    data_file = observe_folder+"bouwens2007_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = (10.**data[:,1])/hubble_h**3.

    data_errorup = (10.**(data[:,1] + data[:,3]) - 10.**data[:,1])/hubble_h**3.
    data_errordown = (10.**data[:,1] - 10.**(data[:,1] + data[:,2]))/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_y/10., fmt='o',label="Bouwens et al. (2007)") 
    return ax

uv_l_z8()


