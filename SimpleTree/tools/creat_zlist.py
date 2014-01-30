import pylab
zlist = pylab.loadtxt("../halofinds")

for z in zlist:
    print 1./(z+1.)
