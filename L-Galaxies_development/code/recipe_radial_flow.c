#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void gas_inflow(int p, double time)
{
	double r_in, r_out, gas_old[RNUM], gasmetal_old[RNUM], newarea, r1, r2, frac, alpha, inflowfrac, rgas, vgas;
	int j, index;
	
	if(Gal[p].ColdGas<1.0e-6) rgas=Gal[p].GasDiskRadius;
   else
  	{
  	  rgas=0.5*radius[0]*Gal[p].ColdGasr[0];
  		for(j=1;j<RNUM;j++) rgas+=(0.5*(radius[j-1]+radius[j])*Gal[p].ColdGasr[j]);
  		rgas=rgas/Gal[p].ColdGas/2.0;      //2.0=mean radius/scale length for exponential disk
    }
   //rgas=Gal[p].GasDiskRadius/3.0;
	 
	//inflowfrac=1.00;
	//alpha=1-1.34e3*time;	
	for (j=0; j<RNUM; j++) 
	{
		inflowfrac=1.00;
		gas_old[j]=Gal[p].ColdGasr[j]*inflowfrac;
		gasmetal_old[j]=Gal[p].MetalsColdGasr[j]*inflowfrac;
		Gal[p].ColdGasr[j]*=(1-inflowfrac);
		Gal[p].MetalsColdGasr[j]*=(1-inflowfrac);
	}
	    	
	for (j=0; j<RNUM; j++) 
	{ 
		vgas=0.7e3; //km/s/Mpc
		alpha=1-vgas*time/Hubble_h; //time unit: (Mpc/h)/(km/s)
		r_out=radius[j]*alpha;
		if(j==0) r_in=0.0;
			else r_in=radius[j-1]*alpha;
		if(r_in<1.0e-8) r_in=1.0e-8;

    /*constant inflow velocity*/			
//		vgas=3.0; // km/s
//		r_out=radius[j]-vgas*time; //time unit: (Mpc/h)/(km/s); radius unit: Mpc/h
//		if(j==0) r_in=0.0;
//	    else r_in=radius[j-1]-vgas*time;
    /*constant inflow velocity*/	
    		
		if(r_out<=radius[0])
			{
				index=0;
				Gal[p].ColdGasr[0]+=gas_old[j];
				Gal[p].MetalsColdGasr[0]+=gasmetal_old[j];
			}
			else
				{
					newarea=r_out*r_out-r_in*r_in;					
					//index=floor(log(1000*r_out)/log(1.2)+3); //The inverse of radius[i]= pow(1.2,i-3)/1000
					for(index=0;radius[index]<=r_out&&index<RNUM;index++);	//ring number of the r_out
					index--;
					r1=radius[index]; r2=r_out;
					while(r1>=r_in)
					{
						frac=(r2*r2-r1*r1)/newarea;
						Gal[p].ColdGasr[index+1]+=frac*gas_old[j];
						Gal[p].MetalsColdGasr[index+1]+=frac*gasmetal_old[j];
						index--;
						r2=r1;
						if(index>-1) r1=radius[index];
							else r1=0.0;				    						
					}
					r1=r_in;
					frac=(r2*r2-r1*r1)/newarea;
					Gal[p].ColdGasr[index+1]+=frac*gas_old[j];
					Gal[p].MetalsColdGasr[index+1]+=frac*gasmetal_old[j];
				}
	}	
}
