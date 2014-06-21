from mass_fn import *
import pylab
import sys
ff = 9
lf = 9
prefix="SA_z6.00"

filter = LGalaxyStruct.properties_used
filter['Sfr'] = True
filter['DiskMass'] = True
filter['BulgeMass'] = True
filter['HotGas'] = True
filter['ColdGas'] = True
filter['BlackholeMass'] = True
filter['Xfrac3d'] = True
file_prefix = "SA_z6.00"
firstfile = 0
lastfile = 127

config = {}
model_names = ["okamoto","noreionization","patchy_I","patcy_II"]
model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization I","Patchy Reionization II"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/","/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/"]
model_plot_patterns = ['r--','b--','g--','m--']

gal = {}
nTrees = {}
nGals = {}
nTreeGals = {}
star = {}
hotgas = {}
coldgas = {}
Blackhole = {}
sfr = {}
for i in range(len(model_names)):
    index = model_names[i]
    (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal(model_paths[i],file_prefix,firstfile,lastfile,filter)
    star[index] = stellar_mass_fn(gal[index],1.,1.e10,50)
    hotgas[index] = hotgas_mass_fn(gal[index],1.e3,1.e12,50)
    coldgas[index] = coldgas_mass_fn(gal[index],1.e3,1.e12,50)
    sfr[index] = sfr_fn(gal[index])



pylab.rc('text', usetex=True)

fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(star[index][0],star[index][1],model_plot_patterns[i],label=model_labels[i])
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Stellar Mass Function z = 6 file "+str(ff)+"-"+str(lf))
#pylab.show()
pylab.savefig('reion_star_'+str(firstfile)+'-'+str(lastfile)+'.pdf',bbox_inches='tight')

# fig = pylab.figure()
# ax = fig.add_subplot(111)
# ax.plot(hot_okamoto_model[0],hot_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
# ax.plot(hot_noreionization_model[0],hot_noreionization_model[1],'b--',label="No Reionization")
# ax.plot(hot_patchyreionization_model[0],hot_patchyreionization_model[1],'m--',label="Patchy Reionization")
# ax.plot(hot_patchyreionization_model_II[0],hot_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
# #ax.plot(hot_patchyreionization_model_III[0],hot_patchyreionization_model_III[1],'k--',label="Patchy Reionization III")
# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)

# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Hot Gass Mass Function z = 6 file "+str(ff)+"-"+str(lf))
# #pylab.show()
# pylab.savefig('reion_hotgas_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


# fig = pylab.figure()
# ax = fig.add_subplot(111)
# ax.plot(cold_okamoto_model[0],cold_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
# ax.plot(cold_noreionization_model[0],cold_noreionization_model[1],'b--',label="No Reionization")
# ax.plot(cold_patchyreionization_model[0],cold_patchyreionization_model[1],'m--',label="Patchy Reionization")
# ax.plot(cold_patchyreionization_model_II[0],cold_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
# #ax.plot(cold_patchyreionization_model_III[0],cold_patchyreionization_model_III[1],'k--',label="Patchy Reionization III")
# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Cold Gass Mass Function z = 6 file "+str(ff)+"-"+str(lf))

# pylab.savefig('reion_coldgas_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


# fig = pylab.figure()
# ax = fig.add_subplot(111)
# ax.plot(bh_okamoto_model[0],bh_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
# ax.plot(bh_noreionization_model[0],bh_noreionization_model[1],'b--',label="No Reionization")
# ax.plot(bh_patchyreionization_model[0],bh_patchyreionization_model[1],'m--',label="Patchy Reionization")
# ax.plot(bh_patchyreionization_model_II[0],bh_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
# #ax.plot(bh_patchyreionization_model_III[0],bh_patchyreionization_model_III[1],'k--',label="Patchy Reionization III")
# #ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Blackhole Mass Function z = 6 file "+str(ff)+"-"+str(lf))

# pylab.savefig('reion_bh_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


# fig = pylab.figure()
# ax = fig.add_subplot(111)
# ax.plot(sfr_okamoto_model[0],sfr_okamoto_model[1],'r--',label="Okamoto et al. (2008)")
# ax.plot(sfr_noreionization_model[0],sfr_noreionization_model[1],'b--',label="No Reionization")
# ax.plot(sfr_patchyreionization_model[0],sfr_patchyreionization_model[1],'m--',label="Patchy Reionization")
# ax.plot(sfr_patchyreionization_model_II[0],sfr_patchyreionization_model_II[1],'g--',label="Patchy Reionization II")
# #ax.plot(sfr_patchyreionization_model_III[0],sfr_patchyreionization_model_III[1],'k--',label="Patchy Reionization III")
# ax.set_yscale("log")
# ax.set_xscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"SFR ($M_\odot/yr$)")
# ax.set_ylabel(r"$N$")
# fig.suptitle("SFR z = 6 file "+str(ff)+"-"+str(lf))

# pylab.savefig('reion_sfr_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')
# #pylab.show()
