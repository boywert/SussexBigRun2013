import pylab
import numpy
import read_lgal
import LGalaxyStruct


def stellar_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['DiskMass']+gadget2msun*gal['BulgeMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(9.0,14.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)

    delta_logM = massftn_x[1]-massftn_x[0]
    pylab.rc('text', usetex=True)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'r-',label="AHF")

    firstfile = 0
    lastfile = 7

     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(subfind_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['DiskMass']+gadget2msun*gal['BulgeMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(9.0,14.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_\star/M_\odot$ $h)$")
    ax.set_ylabel(r"$\phi/$galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_\star/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)


    pylab.show()

stellar_massftn()
