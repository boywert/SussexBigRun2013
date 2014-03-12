
import numpy as np
import matplotlib.pyplot as plt
import pylab
import os
import sys
import math


#import candidates
#sys.path.append('../select_highz/')
sys.path.append('../select_relaxed/')
#sys.path.append('../select_highz_improved')
#import candidates
import rebeccasinourcatalog as candidates
filters=['u','g','r','i1','i2','z','Y','J','H','Ks']




for id in candidates.candidate.keys():
	candidate=candidates.candidate[id]

	upper_limit = False

	m={}
	errm={}
	flux={}
	fluxerr={}
	for f in filters:
		m[f]=candidate[f+'_APER_mag']
		errm[f]=candidate[f+'_APER_magerr']
		
		flux[f]=candidate[f+'_APER_flux']
		fluxerr[f]=candidate[f+'_APER_fluxerr']
		
		snr = flux[f]/fluxerr[f]
		if snr<2.:
			m[f] = -2.5 * np.log10(2.*fluxerr[f]) - 48.6		
			if f=='z': upper_limit = True
		

		
	f1='z'
	f2='Y'
	f3 ='J'
	c1=m[f1]-m[f2]	
	c2=m[f2]-m[f3]
	errc1 = math.sqrt( (errm[f1] ** 2) + (errm[f2] ** 2))
	errc2 = math.sqrt( (errm[f2] ** 2) + (errm[f3] ** 2))
	
	plt.figure()
	plt.ylim([-1.0,10])
	plt.xlim([-1,1])
	plt.xlabel('Y-J')
	plt.ylabel('z-Y')
	plt.axhspan(1.5, 11, xmax=0.75, facecolor='0.5', alpha=0.5)
	
	if upper_limit == True:
		plt.scatter(c2,c1,color='r')
		plt.arrow(c2,c1,0.0,1.0, color='r', fc='r', ec='r', head_width=0.05, head_length=0.2)
	else:
		plt.errorbar(c2,c1,yerr=errc1, xerr=errc2, fmt='ro', marker='.')

	pylab.savefig(id+'.png', bbox_inches=0)
	plt.close()
