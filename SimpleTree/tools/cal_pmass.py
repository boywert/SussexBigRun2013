import pylab

Mpc2m = 3.08567758e22
m2Mpc = 1./Mpc2m
Msun2kg = 1.98855e30
kg2Msun = 1./Msun2kg
m2km = 0.001
h_mass = 1.6737237e-27 #kg
G = 6.674e-11   # SI
H0 = 100.0        # km/s / (Mpc/h)
Msun2Gadget = 1e-10
G = G * (m2km**2.) * (m2Mpc) / (kg2Msun * Msun2Gadget) # (Mpc/h) (km/s)^2 / (1e10Msun/h)
print G
#info
boxsize = 100.0    # Mpc/h
npart =  3072    # per dim
omegam = 0.318  #4320_335
omegab = 0.046
hubble_h = 0.671
pi = pylab.pi

rho_crit_0 = 3.* H0**2 / (8.*pi*G)  # (1e10 Msun/h)/(Mpc/h)^3

pmass = omegam*rho_crit_0*boxsize**3./npart**3.
print "grid mass = ",pmass/8*1e10/hubble_h,"Msun"
print "rho_crit = ",rho_crit_0
print "pmass = ", pmass
print "total H atoms", 1.e10*omegab*rho_crit_0*boxsize**3.*Msun2kg/h_mass/hubble_h
