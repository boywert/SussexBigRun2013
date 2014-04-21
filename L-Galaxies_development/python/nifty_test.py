import matplotlib
matplotlib.use('PDF')
import pylab
import numpy
import read_lgal
import LGalaxyStruct

global subfind_folder
global ahf_folder

subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"


def hotgas_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['HotGas'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['HotGas']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(3.0,15.0))
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
    massf = gadget2msun*gal['HotGas']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(3.0,15.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_{\mathrm{hotgas}}/M_\odot$ $h)$")
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_{\mathrm{hotgas}}/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)
    print "Hot gas mass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM 

    #pylab.show()
    pylab.savefig('hotgas_mass.pdf',bbox_inches='tight')


def coldgas_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['ColdGas'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['ColdGas']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(8.0,13.0))
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
    massf = gadget2msun*gal['ColdGas']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(8.0,13.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_{\mathrm{coldgas}}/M_\odot$ $h)$")
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_{\mathrm{coldgas}}/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)
    print "Cold gas mass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM 


    #pylab.show()
    pylab.savefig('coldgas_mass.pdf',bbox_inches='tight')

def bh_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['BlackHoleMass'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['BlackHoleMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(3.0,10.0))
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
    massf = gadget2msun*gal['BlackHoleMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(3.0,10.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_{\mathrm{bh}}/M_\odot$ $h)$")
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_{\mathrm{bh}}/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    print "Blackhole mass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM 


    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)


    #pylab.show()
    pylab.savefig('bh_mass.pdf',bbox_inches='tight')


def bulge_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['BulgeMass']
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
    massf = gadget2msun*gal['BulgeMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(9.0,14.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_{\mathrm{bulge}}/M_\odot$ $h)$")
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_{\mathrm{bulge}}/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)
    print "Bulge mass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM

    #pylab.show()
    pylab.savefig('bulge_mass.pdf',bbox_inches='tight')



def stellar_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
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
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_\star/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")

    print "Stellar mass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)


    #pylab.show()
    pylab.savefig('stellar_mass.pdf',bbox_inches='tight')

def disk_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = 0
    lastfile = 0

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(ahf_folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['DiskMass']
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
    massf = gadget2msun*gal['DiskMass']
    mass = numpy.log10(massf)
    stellarmass = pylab.histogram(mass,bins=20,range=(9.0,14.0))
    print stellarmass
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    ax.set_xlabel(r"$\log(M_{\mathrm{disk}}/M_\odot$ $h)$")
    ax.set_ylabel(r"galaxies$/(Mpc^3 h^{-3})/\Delta \log(M_{\mathrm{disk}}/M_\odot$ $h)$")
    ax.plot(massftn_x,massftn_y/boxsize**3./delta_logM,'b-',label="SUBFIND")
    print "Diskmass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM/boxsize**3./delta_logM
    ax.set_yscale("log")
    ax.legend(loc='upper right',ncol=1, fancybox=True)


   # pylab.show()                
    pylab.savefig('disk_mass.pdf',bbox_inches='tight')

stellar_massftn()
bulge_massftn()
disk_massftn()
bh_massftn()
hotgas_massftn()
coldgas_massftn()
