import numpy
import read_lgal
import LGalaxyStruct

def stellar_mass_fn(folder,file_prefix,firstfile,lastfile,mass_min=1.,mass_max=1.e20,nbins=20):
    print "reading "+folder
    gadget2msun=10.e10
    boxsize = 47.0
    hubble_h = 0.7
    filter = LGalaxyStruct.properties_used
    filter['DiskMass'] = True
    filter['BulgeMass'] = True

    (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal(folder,file_prefix,firstfile,lastfile,filter)
    massf = gadget2msun*gal['DiskMass']+gadget2msun*gal['BulgeMass']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)


