from stellar_mass_fn import *
import pylab
import sys
ff = 3
lf = 127
prefix="SA_z6.00"
okamoto_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/",prefix,ff,lf,1.e3,1.e12,50)
noreionization_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/",prefix,ff,lf,1.e3,1.e12,50)
#fullreionization_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/full_reionization/","SA_z6.00",ff,lf,1.e3,1.e12,50)
patchyreionization_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization/",prefix,ff,lf,1.e3,1.e12,50)
patchyreionization_model_II = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/",prefix,ff,lf,1.e3,1.e12,50)
print patchyreionization_model
pylab.rc('text', usetex=True)
fig = pylab.figure()
ax = fig.add_subplot(111)


ax.plot(okamoto_model[0],okamoto_model[1],'r-',label="Okamoto et al. (2008)")
ax.plot(noreionization_model[0],noreionization_model[1],'b--',label="No Reionization")
#ax.plot(fullreionization_model[0],fullreionization_model[1],'g-',label="Full Reionization")
ax.plot(patchyreionization_model[0],patchyreionization_model[1],'m--',label="Patchy Reionization")
ax.plot(patchyreionization_model_II[0],patchyreionization_model_II[1],'g--',label="Patchy Reionization II")

ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Stellar Mass Function z = 6")
pylab.savefig('reion_model_8.pdf',bbox_inches='tight')
