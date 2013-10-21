import numpy
import os
import math

global SNAPfile
global AHFdir
global AHFprefix
global G
global m2Mpc
global m2km
global kpc2Mpc
global Msun2Gadget
global kg2Msun

G = 6.67384e-11 # m^3/(kgs^2)
m2Mpc = 1./3.08567758e22
m2km = 0.001
kpc2Mpc = 0.001
Msun2Gadget = 1.e-10
kg2Msun = 1.989e-30
#change to (Mpc/h) (km/s)^2 / (1e10Msun/h)
G = G*m2Mpc*m2km**2./(Msun2Gadget*kg2Msun)

AHFdir = "/scratch/datasetI"
AHFprefix = "62.5_dm"
SUSSINGtree = "/export/research/virgo/Boyd/SUSSING2013/DATASET_I/MergerTree"
SNAPfile = "/scratch/datasetI/data_snaplist.txt"


def readAHFascii():
    halocat = {}
    timesnap = numpy.loadtxt(SNAPfile)
    for time in timesnap:
        zstring = "%.3f" % (time[2])
        #print zstring[len(zstring)-1]
        filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        #print "checking "+filename
        if os.path.isfile(filename) == False:
            zstring = "%.3f" % (time[2]+0.00001)
            #print zstring[len(zstring)-1] 
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            zstring = "%.3f" % (time[2]-0.00001)
            #print "checking "+filename
            #print zstring[len(zstring)-1] 
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            print "checking error"+filename
            exit()
        print "Reading "+filename
        stat = os.stat(filename)
        #print stat.st_size
        if(stat.st_size > 384):
            data = numpy.loadtxt(filename)
            for halo in data:
                hid = long(halo[0])
                halocat[hid] = {}
                halocat[hid]["ID"] = hid
                halocat[hid]["Mvir"] = halo[3]*Msun2Gadget
                halocat[hid]["Len"] = halo[4]
                halocat[hid]["Pos"] = (halo[5]*kpc2Mpc,halo[6]*kpc2Mpc,halo[7]*kpc2Mpc)
                halocat[hid]["Vel"] = (halo[8],halo[9],halo[10])
                halocat[hid]["Vmax"] = halo[16]
                halocat[hid]["VelDisp"] = halo[18]
                # use Peebles lambdaE definition to find angular momentum 
                halocat[hid]["Ep"] = halo[38]
                halocat[hid]["Ek"] = halo[39]
                total_energy = math.fabs((halo[38] + halo[39])*Msun2Gadget)
                halocat[hid]["LambdaE"] = halo[20]
                J = halo[20]*G*halocat[hid]["Mvir"]**(3./2.)/total_energy**(0.5)
                halocat[hid]["TotalEnergy"] = total_energy
                halocat[hid]["Spin"] = (halo[21]*J,halo[22]*J,halo[23]*J)
                halocat[hid]["FirstProgenitor"] = -1
                halocat[hid]["NextProgenitor"] = -1
                halocat[hid]["Descendant"] = -1
                halocat[hid]["HostHalo"] = long(halo[1])
                if(halocat[hid]["HostHalo"] == 0):
                    halocat[hid]["HostHalo"] = -1
                halocat[hid]["NextHalo"] = -1
                halocat[hid]["SnapNum"] = long(time[0])
    return halocat

def makeStuctree(halocat):
    timesnap = numpy.loadtxt(SNAPfile)
    for haloc in halocat:
        halo = halocat[haloc]
        hosthalo= halo["HostHalo"]
        if(hosthalo > -1):
            cursub = hosthalo
            while cursub > -1:
                curid = halocat[cursub]["ID"]
                cursub = halocat[cursub]["NextHalo"]
            print "change",halocat[curid]["NextHalo"],"to",halloc
            halocat[curid]["NextHalo"] = haloc
        
    return halocat

def readSussingtree(SUSSINGtree,halocat):
    f = open(SUSSINGtree)
    line = f.read().splitlines()
    count = 0;
    for (i,item) in enumerate(line):
        if(i == 2):
            totalhalo = long(item)
            print "tree",totalhalo,"halocat",len(halocat)
        if(i >= 3):
            col = item.split()
            if(count == 0):
                if(col[0] == "END"):
                    return
                if(len(col) != 2):
                    print "line",i,"has error"
                    exit()
                else:
                    haloid = col[0]
                    nprog = long(col[1])
                    count = nprog
            else:
                if(len(col) != 1):
                    print "line",i,"has error"
                    exit()
                else:
                    progid = long(col[0])
                    count -= 1

    return halocat
