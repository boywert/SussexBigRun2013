;=========================================================================
  
   PRO plot_lf2,mag,mag2,magname=magname,$
               ps=ps,box=box,Hubble_h=Hubble_h,$
               min=min,max=max,binsize=binsize

;=========================================================================
;
;  Example procedure to create luminosity function from magnitudes
;
;-------------------------------------------------------------------------

; Inputs:
; mag=array of magnitudes

magname = "M!DUV=GFUV!n"
ps = ps

;; Dust
Ngals=n_elements(mag)
if Ngals eq 0 then message,'plot_lf called with zero-length magnitude array'
if n_elements(magname) eq 0 then magname='magnitude'
if n_elements(box) eq 0 then begin
   print,'Using boxsize=500 Mpc'
   box=500.
endif
if n_elements(Hubble_h) eq 0 then begin
   print,'Using Hubble_h= 0.704'
   Hubble_h= 0.704
endif
if n_elements(min) eq 0 then min = -30.
if n_elements(max) eq 0 then max = 0.
if n_elements(binsize) eq 0 then binsize = 0.25

; volume in Mpc
volume = (box/Hubble_h)^3.0
log_hubble = alog10(Hubble_h)


;-------------------------------------------------------------------------------

;if (keyword_set(ps)) then begin
    set_plot, 'PS'
    device, filename = 'lf.ps', xsize = 12, ysize = 11, /color, $
      xoffset=1, $
      yoffset=5
;endif
  
; Create histogram of data and data points
  hist=(histogram(mag,locations=c,min=min,max=max,binsize=binsize))
  ; x is centre of bins
  x=c+0.5*binsize
  ; y is normalised to volume and binsize
  y=hist/(volume*binsize)

  ;without dust
  hist=(histogram(mag2,locations=c,min=min,max=max,binsize=binsize))
  ; x is centre of bins
  x2=c+0.5*binsize
  ; y is normalised to volume and binsize
  y2=hist/(volume*binsize)


; h
  x = x - (5*ALOG10(Hubble_h))
  y = y/(Hubble_h^3)
  x2 = x2 - (5*ALOG10(Hubble_h))
  y2 = y2/(Hubble_h^3)

; Set axes
  xmin=min-0.5
  xmax=max+0.5
  index=where(y gt 1e-30)
  ylogmin=-10 > floor(min(alog10(y[index])))
  ylogmax=ceil(max(alog10(y[index])))

; New axis limits
  xmin = -23
  xmax = -16
  ylogmin = -6
  ylogmax = -1
  
; Legend
;  meanings = ['Straight','Dashed']
;  psyms=[-1,-2,-4]
;  lines=[1,2,4]

;Obs data
;z=4
;z=5
  b07z5x = [-22.66177,-22.17646,-21.66177,-21.17647,-20.66177, -20.17647, -19.67647, -19.17647, -18.66177, -18.16177, -17.66177, -17.16177]
  b07z5y = [-5.32755, -5.32755, -4.46068, -3.94056, -3.50712, -2.94984, -2.75170, -2.68978, -2.29350, -2.28111, -2.10774]
;z=6
  b07z6x = [-22.13235, -21.63235, -21.14706, -20.63235, -20.13235, -19.63235, -18.89706, -17.88235]
  b07z6y = [-5.16656, -5.29040, -4.15108, -3.89102, -3.28421, -3.09845, -2.71455, -2.26873]

;z=7
  b11x = [-20.80656,-20.30451,-19.81271,-19.31066,-18.80861,-18.30656,-17.80451]
  b11y = [-3.78723,-3.84043,-3.29787,-2.97872,-2.62766,-2.46808,-2.17021]

  mc10x = [-20.50062, -19.75113, -19.25832, -18.65257, -18.14949]
  mc10y = [-3.79200, -3.30133, -3.10933, -2.96000, -2.70400]

  oe10x = [-18.35482, -19.15565, -19.95647, -20.75729]
  oe10y = [-2.57600, -2.91733, -3.33333, -4.00544]
;z=8
  b10z8x=[-18.87013, -18.37662]
  b10z8y=[-2.95372, -3.05289]

  b11z8x =[-20.745,-20.147, -19.549, -18.938, -18.34, -17.742]
  b11z8y =[-3.95537, -3.60826, -3.4, -3.00331, -2.81488, -2.34876]

  mc10z8x = [-19.79221, -18.79221]
  mc10z8y = [-3.74711, -3.36033]


; Plot
	redshift="z10.94"

  plot,x2, y2, /ylog, title = 'Hen12 uv      '+redshift, $
        xstyle=1, xrange = [xmin,xmax], xtitle = magname,$
        ystyle=1, yrange = 10.^[ylogmin,ylogmax],$
        ytitle = '!4U!3 /(Mpc!u-3!nmag!u-1!n)'
  oplot, x,y, linestyle = 2, color = 2
  oplot, b07z5x-5*ALOG10(hubble_h), (10^(b07z5y))/hubble_h^3, psym=2, color=4
 ; oplot, b07z6x, 10^(b07z6y), psym=2, color=4
 ; oplot, b11x, 10^(b11y), color = 4, psym = 2
 ; oplot, mc10x, 10^(mc10y), psym = 7, color=4
 ; oplot, oe10x, 10^(oe10y), psym = 1, color=4
 ; oplot, b10z8x, 10^(b10z8y) ,psym=2,color=4
 ; oplot, b11z8x, 10^(b11z8y),psym=5, color=4
 ; oplot, mc10z8x, 10^(mc10z8y),psym=4, color=4

  ;xyouts, 0.5, 0.5,/normal, 'Sine Curve'
 ; xyouts, -22.5, 0.02, 'No Dust'
 ; plots, [-21.25,-20.5], [0.02, 0.02], linestyle=0

 ; xyouts, -22.5, 0.01, 'With Dust'
 ; plots, [-21.25,-20.5], [0.01, 0.01], linestyle=2, color=2;

 ; xyouts, -22.5, 0.005, 'Bouwens 2007'
 ; plots, -20.15, 0.005, psym=2, color=4

 ; xyouts, -22.5, 0.002, 'Oesch 2010'
 ; plots, -20.15, 0.002, psym=1, color=4

 ; xyouts, -22.5, 0.001, 'Bouwens 2011'
 ; plots, -20.15, 0.001, psym=2, color=4
 

;print to file

;openw,1,redshift+'_sfr_x.cat' ;open the file to write
;printf,1,"# 1 x"
;printf,1,x
;close, 1 ;close file
;
;openw,2,redshift+'_sfr_y.cat'
;printf,2,"# 1 y"
;printf,2,y
;close,2

;openw,3,redshift+'_sfr_x.cat' ;open the file to write
;printf,3,"# 1 x"
;printf,3,x
;close, 3 ;close file
;
;openw,4,redshift+'_sfr_y.cat'
;printf,4,"# 1 y"
;printf,4,y
;close,4

  

;if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
;endif

;------------------------------------------------------------------------------


end
