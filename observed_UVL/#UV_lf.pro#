pro plot_UV_LF_highz, Datadir, hubble_h, volume_MR, volume_MRII, $
                      G0_MR, G1_MR, G2_MR, G3_MR, G0_MRII, G1_MRII, G2_MRII, G3_MRII, $
                      StellarMass_MR_0, StellarMass_MR_1, StellarMass_MR_2, StellarMass_MR_3, $
                      StellarMass_MRII_0, StellarMass_MRII_1, StellarMass_MRII_2, StellarMass_MRII_3

;NIC_F160W


  bin=0.5
  xmin=-23.0
  xmax=-16.0
  ymin=1.e-6
  ymax=1.e-1
 
for i=0, 3 do begin

   if(i eq 0) then begin
      G=G0_MRII 
      plot_obs_uv_z5, Datadir, xmin, xmax, ymin, ymax
   endif
   
   if(i eq 1) then begin
      G=G1_MRII
      multiplot
      plot_obs_uv_z6, Datadir, xmin, xmax, ymin, ymax
   endif

   if(i eq 2) then begin
      G=G2_MRII 
      multiplot
      plot_obs_uv_z7, Datadir, xmin, xmax, ymin, ymax
   endif
   
   if(i eq 3) then begin
      G=G3_MRII
      multiplot
      plot_obs_uv_z8, Datadir, xmin, xmax, ymin, ymax
   endif

;G11
  ;readcol,Datadir+'guo_wmp7_bband_m05_z10.txt',mag,phi
  ;oplot,mag,phi, color=70,linestyle=2, thick=3

;NEW MODEL
  mag=G.MagDust[31]
  ;mag=G.MagDust[0]
  hist=histogram(mag,locations=c,min=xmin-5.0,max=xmax+5.0,binsize=bin)
  oplot,c+bin/2.0-5*Alog10(hubble_h),hist/(volume_MRII*bin), color=70, thick=6
   

endfor

end

pro plot_obs_uv_z5, Datadir, xmin, xmax, ymin, ymax

  hubble_h=0.7

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, /ylog, ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'

  readcol,Datadir+'bouwens2007_z5.txt',mag,phi,errordown,errorup
  symbols, 2, 0.7  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  
  
  oploterror, [-22.5,-22.5], [10^(-1.28),10^(-1.28)], [0.012,0.012], $
              psym=8, color=200, errcolor=200, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.35), 'Bouwens (2007)'

  xyouts, -18., 1.0e-4, 'z=5',charsize=1.3,charthick=4

end




pro plot_obs_uv_z6, Datadir, xmin, xmax, ymin, ymax

  hubble_h=0.7

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, /ylog      

  readcol,Datadir+'bouwens2007_z6.txt',mag,phi,errordown,errorup
  symbols, 2, 0.7
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  
  
  oploterror, [-22.5,-22.5], [10^(-1.28),10^(-1.28)], [0.012,0.012], $
              psym=8, color=200, errcolor=200, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.35), 'Bouwens (2007)'
  
  xyouts, -18., 1.0e-4, 'z=6',charsize=1.3,charthick=4

end




pro plot_obs_uv_z7, Datadir, xmin, xmax, ymin, ymax

  hubble_h=0.7

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, /ylog, $
        xtitle = 'M!D1600!N(AB)-5log!D10!Nh', ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'

  readcol,Datadir+'bouwens2011_z7.txt',mag,phi,eleft,eright,errordown,errorup
  symbols, 2, 0.7  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  
  
  readcol,Datadir+'mclure2010_z7.txt',mag,phi,errordown,errorup
  symbols, 30, 0.3  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 120, $
              errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 120, $
              errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
  
  readcol,Datadir+'oesch2010_z7.txt',mag,phi,errordown,errorup
  symbols, 32, 0.25  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 0, $
              errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 0, $
              errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
   
  
  symbols, 2, 0.7
  oploterror, [-22.5,-22.5], [10^(-1.28),10^(-1.28)], [0.012,0.012], $
              psym=8, color=200, errcolor=200, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.35), 'Bouwens (2011)'
  
  
  symbols, 30, 0.3
  oploterror, [-22.5,-22.5], [10^(-1.55),10^(-1.55)], [0.005,0.005], $
              psym=8, color=120, errcolor=120, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.6), 'McLure (2010)'
  
  symbols, 32, 0.25 
  oploterror, [-22.5,-22.5], [10^(-1.8),10^(-1.8)], [0.003,0.003], $
              psym=8, color=0, errcolor=0, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.85), 'Oesch (2010)'

  xyouts, -18., 1.0e-4, 'z=7',charsize=1.3,charthick=4

end




pro plot_obs_uv_z8, Datadir, xmin, xmax, ymin, ymax

  hubble_h=0.7

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, /ylog, $
        xtitle = 'M!D1600!N(AB)-5log!D10!Nh'

  readcol,Datadir+'bouwens2011_z8.txt',mag,phi,eleft,eright,errordown,errorup
  symbols, 2, 0.7  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 200, $
              errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  
  
  readcol,Datadir+'mclure2010_z8.txt',mag,phi,errordown,errorup
  symbols, 30, 0.3  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 120, $
              errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 120, $
              errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
  
  readcol,Datadir+'bouwens2010_z8.txt',mag,phi,errordown,errorup
  symbols, 32, 0.25  
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^(phi+errorup)-10^phi)/hubble_h^3, color = 0, $
              errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(hubble_h), (10^phi)/hubble_h^3, (10^phi-10^(phi+errordown))/hubble_h^3, color = 0, $
              errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
  
  
  
  oploterror, [-22.5,-22.5], [10^(-1.28),10^(-1.28)], [0.012,0.012], $
              psym=8, color=200, errcolor=200, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.35), 'Bouwens (2011)'
  
  
  
  symbols, 30, 0.3
  oploterror, [-22.5,-22.5], [10^(-1.55),10^(-1.55)], [0.005,0.005], $
              psym=8, color=120, errcolor=120, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.6), 'McLure (2010)'
  
  symbols, 32, 0.25 
  oploterror, [-22.5,-22.5], [10^(-1.8),10^(-1.8)], [0.003,0.003], $
              psym=8, color=0, errcolor=0, HATLENGTH = 80.0
  xyouts, -22.25, 10^(-1.85), 'Bouwens (2010)'


  xyouts, -18., 1.0e-4, 'z=8',charsize=1.3,charthick=4

end
