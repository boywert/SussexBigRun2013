import numpy
import os

AHFdir = "/scratch/datasetI"
AHFprefix = "62.5_dm"
SUSSINGtree = "/export/research/virgo/Boyd/SUSSING2013/DATASET_I/MergerTree"
SNAPfile = "/scratch/datasetI/data_snaplist.txt"
halocat = {}
global halocat
def readAHFascii(SNAPfile,AHFdir,AHFprefix):
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
                halocat[halo[0]] = {}
                halocat[halo[0]]["Mvir"] = halo[3]

def readSussingtree(SUSSINGtree):
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
                    count--

                
            
readAHFascii(SNAPfile,AHFdir,AHFprefix)
readSussingtree(SUSSINGtree)
