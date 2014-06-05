from stellar_mass_fn import *
import pylab

okamoto_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/","SA_z6.00",0,127,1.e3,1.e12,50)

print okamoto_model
pylab.rc('text', usetex=True)
fig = pylab.figure()
ax = fig.add_subplot(111)


ax.plot(okamoto_model[0],okamoto_model[1],'r-',label="Okamoto et al. (2008)")



ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log{M/M_\odot h}$")
ax.set_ylabel("STMF")
fig.suptitle("Stellar Mass Function z = 6")
pylab.show()
