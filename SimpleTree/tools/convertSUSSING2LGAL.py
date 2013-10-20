import numpy
import os

AHFdir = "/scratch/datasetI"
AHFprefix = "62.5_dm"
SUSSINGtree = ""
SNAPfile = "/scratch/datasetI/data_snaplist.txt"

def readAHFascii(SNAPfile,AHFdir,AHFprefix):
    timesnap = numpy.loadtxt(SNAPfile)
    halocat = {}
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
                halocat[halo[0]] = {}
                halocat[halo[0]]["Mvir"] = halo[3]
                halocat[halo[0]]["Len"] = int(halo[2])
readAHFascii(SNAPfile,AHFdir,AHFprefix)
