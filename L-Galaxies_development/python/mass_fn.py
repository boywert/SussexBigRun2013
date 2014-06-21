import numpy
import read_lgal
import LGalaxyStruct
gadget2msun=10.e10
boxsize = 47.0
hubble_h = 0.7

def stellar_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = gadget2msun*gal['DiskMass']+gadget2msun*gal['BulgeMass']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def hotgas_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = gadget2msun*gal['HotGas']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def coldgas_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = gadget2msun*gal['ColdGas']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def bh_mass_fn(gal,mass_min=1.0,mass_max=1.e8,nbins=50):
    massf = numpy.log10(gal['BlackHoleMass'])
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    for i in range(len(massftn_x)):
        massftn_x[i] = pow(10.,massftn_x[i])
    return (massftn_x,massftn_y)


def sfr_fn(gal,mass_min=0.000001,mass_max=1.,nbins=50):
    massf = numpy.log10(gal['Sfr'])
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    for i in range(len(massftn_x)):
        massftn_x[i] = pow(10.,massftn_x[i])
    return (massftn_x,massftn_y)
