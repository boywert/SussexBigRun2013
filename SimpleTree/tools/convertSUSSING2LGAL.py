import numpy
import os

AHFdir = "/export/data/virgo/SUSSING2013_DATA/datasetI"
AHFprefix = "62.5_dm"
SUSSINGtree = ""
SNAPfile = "/export/data/virgo/SUSSING2013_DATA/datasetI/data_snaplist.txt"

def readAHFascii(SNAPfile,AHFdir,AHFprefix):
    timesnap = numpy.loadtxt(SNAPfile)
    halocat = {}
    for time in timesnap:
        filename = "%s/%s_%03d.z%3.3f.AHF_halos" % (AHFdir, AHFprefix, time[0], time[2])
        print "Reading "+filename
        stat = os.stat(filename)
        print stat.st_size
        #data = numpy.loadtxt(filename)
        #for halo in data:
        #    print halo

readAHFascii(SNAPfile,AHFdir,AHFprefix)
