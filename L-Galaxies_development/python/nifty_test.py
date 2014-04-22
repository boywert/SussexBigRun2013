import matplotlib
matplotlib.use('PDF')
import pylab
import numpy
import read_lgal
import LGalaxyStruct

global folder1
global folder2
global first1
global last1
global fitst2
global last2
global file_prefix
global boxsize1
global boxsize2
global label1
global label2

folder1 = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
folder2 = "/mnt/lustre/scratch/cs390/nIFTy/500_test/outputs/"
first1  = 0
last1 = 7
first2 = 10
last2 = 10
file_prefix = "SA_z0.00"
boxsize1 = 62.5
boxsize2 = 500.0/8.0
label1 = "62.5Mpc/h SF"
label2 = "500Mpc/h SF"

def hotgas_massftn():
    gadget2msun=10.e10
    boxsize = 47.0
    max_mag=-16.0
    min_mag = -23.
    nbins=14
    hubble_h = 0.7
    #subfind_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dmSF/outputs/"
    #ahf_folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/outputs/"
    
    firstfile = first1
    lastfile = last1

    filter = LGalaxyStruct.properties_used
    filter['HotGas'] = True

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

    firstfile = first2
    lastfile = last2

     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)
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
    
    firstfile = first1
    lastfile = last1

    filter = LGalaxyStruct.properties_used
    filter['ColdGas'] = True

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

 
     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)
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
    
    firstfile = first1
    lastfile = last1

    filter = LGalaxyStruct.properties_used
    filter['BlackHoleMass'] = True

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

    firstfile = first2
    lastfile = last2
     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)
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
    
    firstfile = first1
    lastfile = last1 

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

    firstfile = first2
    lastfile = last2

     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)
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
    
    firstfile = first1
    lastfile = last1

    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

    firstfile = first2
    lastfile = last2

     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)

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

    #file_prefix = "SA_z0.00"    
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder1,file_prefix,first1,last1,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize1**3./delta_logM,'r-',label=label1)

    firstfile = 0
    lastfile = 7

     
    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder2,file_prefix,first2,last2,filter)
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
    ax.plot(massftn_x,massftn_y/boxsize2**3./delta_logM,'b-',label=label2)
    print "Diskmass"
    for i in range(len(massftn_x)):
        print massftn_x[i],"\t",massftn_y[i]/boxsize**3./delta_logM
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
