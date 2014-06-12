from mass_fn import *
import pylab
import sys
ff = 10
lf = 10
prefix="SA_z6.00"
star_okamoto_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/",prefix,ff,lf,1.e3,1.e12,50)
star_noreionization_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/",prefix,ff,lf,1.e3,1.e12,50)
star_patchyreionization_model = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/",prefix,ff,lf,1.e3,1.e12,50)
star_patchyreionization_model_II = stellar_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/",prefix,ff,lf,1.e3,1.e12,50)


hot_okamoto_model = hotgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/",prefix,ff,lf,1.e3,1.e12,50)
hot_noreionization_model = hotgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/",prefix,ff,lf,1.e3,1.e12,50)
hot_patchyreionization_model = hotgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/",prefix,ff,lf,1.e3,1.e12,50)
hot_patchyreionization_model_II = hotgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/",prefix,ff,lf,1.e3,1.e12,50)


cold_okamoto_model = coldgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/",prefix,ff,lf,1.e3,1.e12,50)
cold_noreionization_model = coldgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/",prefix,ff,lf,1.e3,1.e12,50)
cold_patchyreionization_model = coldgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/",prefix,ff,lf,1.e3,1.e12,50)
cold_patchyreionization_model_II = coldgas_mass_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/",prefix,ff,lf,1.e3,1.e12,50)

sfr_okamoto_model = sfr_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/",prefix,ff,lf)
sfr_noreionization_model = sfr_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/",prefix,ff,lf)
sfr_patchyreionization_model = sfr_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/",prefix,ff,lf)
sfr_patchyreionization_model_II = sfr_fn("/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/",prefix,ff,lf)

pylab.rc('text', usetex=True)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(star_okamoto_model[0],star_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
ax.plot(star_noreionization_model[0],star_noreionization_model[1],'b--',label="No Reionization")
ax.plot(star_patchyreionization_model[0],star_patchyreionization_model[1],'m--',label="Patchy Reionization")
ax.plot(star_patchyreionization_model_II[0],star_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Stellar Mass Function z = 6 file "+str(ff)+"-"+str(lf))
#pylab.show()
pylab.savefig('reion_star_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(hot_okamoto_model[0],hot_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
ax.plot(hot_noreionization_model[0],hot_noreionization_model[1],'b--',label="No Reionization")
ax.plot(hot_patchyreionization_model[0],hot_patchyreionization_model[1],'m--',label="Patchy Reionization")
ax.plot(hot_patchyreionization_model_II[0],hot_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Hot Gass Mass Function z = 6 file "+str(ff)+"-"+str(lf))
#pylab.show()
pylab.savefig('reion_hotgas_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(cold_okamoto_model[0],cold_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
ax.plot(cold_noreionization_model[0],cold_noreionization_model[1],'b--',label="No Reionization")
ax.plot(cold_patchyreionization_model[0],cold_patchyreionization_model[1],'m--',label="Patchy Reionization")
ax.plot(cold_patchyreionization_model_II[0],cold_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Cold Gass Mass Function z = 6 file "+str(ff)+"-"+str(lf))

pylab.savefig('reion_coldgas_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(sfr_okamoto_model[0],sfr_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
ax.plot(sfr_noreionization_model[0],sfr_noreionization_model[1],'b--',label="No Reionization")
ax.plot(sfr_patchyreionization_model[0],sfr_patchyreionization_model[1],'m--',label="Patchy Reionization")
ax.plot(sfr_patchyreionization_model_II[0],sfr_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
ax.set_yscale("log")
ax.set_xscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
ax.set_xlabel(r"SFR")
ax.set_ylabel(r"$N$")
fig.suptitle("SFR z = 6 file "+str(ff)+"-"+str(lf))

pylab.savefig('reion_sfr_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')
#pylab.show()
