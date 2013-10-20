import numpy

AHFdir = "/export/data/virgo/SUSSING2013_DATA/datasetI/"
AHFprefix = "62.5_dm_"
SUSSINGtree = ""
SNAPfile = "/export/data/virgo/SUSSING2013_DATA/datasetI/data_snaplist.txt"

def readAHFascii(SNAPfile,AHFdir,AHFprefix):
    timesnap = numpy.loadtxt(SNAPfile)
    for time in timesnap:
        print time
