from mass_fn import *
import pylab
import sys
import matplotlib.pyplot as plt

filter = LGalaxyStruct.properties_used
filter['Sfr'] = True
filter['DiskMass'] = True
filter['BulgeMass'] = True
filter['HotGas'] = True
filter['ColdGas'] = True
filter['BlackholeMass'] = True
filter['Xfrac3d'] = True
filter['HaloM_Crit200'] = True
file_prefix = "SA_z8.06"
firstfile = 1
lastfile = 10

config = {}
model_names = ["okamoto","noreionization","patchy_I","patchy_II"]
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
fig.suptitle("Stellar Mass Function z = 6 file "+str(firstfile)+"-"+str(lastfile))
pylab.savefig('reion_star_'+str(firstfile)+'-'+str(lastfile)+'.pdf',bbox_inches='tight')


fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(hotgas[index][0],hotgas[index][1],model_plot_patterns[i],label=model_labels[i])

ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Hot Gass Mass Function z = 6 file "+str(firstfile)+"-"+str(lastfile))
pylab.savefig('reion_hotgas_'+str(firstfile)+'-'+str(lastfile)+'.pdf',bbox_inches='tight')


fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(coldgas[index][0],coldgas[index][1],model_plot_patterns[i],label=model_labels[i])

ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Cold Gass Mass Function z = 6 file "+str(firstfile)+"-"+str(lastfile))

pylab.savefig('reion_coldgas_'+str(firstfile)+'-'+str(lastfile)+'.pdf',bbox_inches='tight')


# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(Blackhole[index][0],Blackhole[index][1],model_plot_patterns[i],label=model_labels[i])

# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Blackhole Mass Function z = 6 file "+str(ff)+"-"+str(lf))

# pylab.savefig('reion_bh_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


fig = plt.figure()
ax = fig.add_subplot(111)

z = gal["okamoto"]["Xfrac3d"]
x = numpy.log10(gal["okamoto"]["HaloM_Crit200"]*1.e10)
y = numpy.log10((gal["okamoto"]["BulgeMass"]+gal["okamoto"]["DiskMass"])*1.e10)

ax.scatter(x, y, c=z, s=2, cmap=pylab.cm.cool, edgecolors='None')
plt.colorbar()
# ax.scatter(x,y,s=2)
# ax.scatter(range(len(gal["okamoto"]["DiskMass"])),(gal["okamoto"]["BulgeMass"]+gal["okamoto"]["DiskMass"])*10.e10,s=2,color='red',label='oka')
# ax.scatter(range(len(gal["patchy_II"]["DiskMass"])),(gal["patchy_II"]["BulgeMass"]+gal["patchy_II"]["DiskMass"])*10.e10,s=2,color='green',label='patchy')

#ax.set_yscale("log")
f = open("test","w+")
for i in range(len(gal["okamoto"]["Xfrac3d"])):
    f.write("%f %f\n" % (gal["okamoto"]["Xfrac3d"][i],(gal["okamoto"]["BulgeMass"][i]+gal["okamoto"]["DiskMass"][i])*10.e10))
f.close()
#leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# #leg.get_frame().set_linewidth(0)
# ax.set_ylabel(r"SFR ($M_\odot/yr$)")
# ax.set_xlabel(r"$\log(M/M_\odot h)$")

# fig.suptitle("SFR z = 6 file "+str(firstfile)+"-"+str(lastfile))
#fig.show()
pylab.savefig('reion_sfr_'+str(firstfile)+'-'+str(lastfile)+'.png',bbox_inches='tight')

