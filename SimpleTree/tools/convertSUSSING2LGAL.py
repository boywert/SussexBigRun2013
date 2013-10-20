import numpy
import os

AHFdir = "/export/data/virgo/SUSSING2013_DATA/datasetI"
AHFprefix = "62.5_dm"
SUSSINGtree = ""
SNAPfile = "/export/data/virgo/SUSSING2013_DATA/datasetI/data_snaplist.txt"
def properround(number,digits):
    strfmt = "%."+str(digits+1)+"f"
    strnum = strfmt % (number)
    if strnum[len(strnum)-1] == "5":
        if long(strnum[len(strnum)-2])%2 == 0:
            strout = strnum[0:(len(strnum)-1)]
        else:
            strnum[len(strnum)-2] = str(long(strnum[len(strnum)-2])+1)
            strout = strnum[0:(len(strnum)-1)]
    else:
        strfmt = "%."+str(digits)+"f"
        strout = strfmt % (number)
    return strout
def readAHFascii(SNAPfile,AHFdir,AHFprefix):
    timesnap = numpy.loadtxt(SNAPfile)
    halocat = {}
    for time in timesnap:
        zstring = "%.3f" % (time[2])
        filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        print "checking "+filename
        if os.path.isfile(filename) == False:
            print "checking "+filename
            zstring = "%.3f" % (time[2]+0.001)
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            print "checking "+filename
            zstring = "%.3f" % (time[2]-0.001)
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            print "checking "+filename
            exit()
        print "Reading "+filename
        stat = os.stat(filename)
        print stat.st_size
        #data = numpy.loadtxt(filename)
        #for halo in data:
        #    print halo

readAHFascii(SNAPfile,AHFdir,AHFprefix)
