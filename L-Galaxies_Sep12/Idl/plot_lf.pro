;=========================================================================
  
   PRO plot_lf,mag,magname=magname,$
               ps=ps,box=box,Hubble_h=Hubble_h,LFbinwifth=LFbinwidth

;=========================================================================
;
;  Example procedure to create historgram (mass funciton or luminosity function)
;
;-------------------------------------------------------------------------

; Inputs:
; mag=array of magnitudes
if n_elements(magname) eq 0 then magname='magnitude'
if n_elements(box) eq 0 then box=500.
if n_elements(Hubble_h) eq 0 then Hubble_h= 0.704
if n_elements(LFbinwidth) eq 0 then LFbinwidth = 0.25

; volume in Mpc
volume = (box/Hubble_h)^3.0
log_hubble = alog10(Hubble_h)

;-------------------------------------------------------------------------------

; An example of how to create data for a plot, in this case
; luminosities

Ngals=n_elements(Type)

; Calculate luminosities
L_B  = fltarr(Ngals)
L_K  = fltarr(Ngals)
L_bj = fltarr(Ngals)
Mag_bj = fltarr(Ngals)
B_cut = 50.0
V_cut = 50.0
K_cut = 50.0
bj_cut = -16.0                  ; bj magnitude limit in survey
mbj_sun = 5.33                  ; magnitude of sun in bj band
mk_sun = 3.3                    ; magnitude of sun in K band
mb_sun = 5.48                   ; magnitude of sun in B band
L_bjmin = 10.0^((mbj_sun - bj_cut) / 2.5)
L_Bmin = 1.e6
L_Kmin = 1.e6
type_pivot = 0.8                ; type_pivot<0.8 ---> blue
                                ; type_pivot>0.8 ---> red

; Calculate luminosities and bj mag
for i=0L, Ngals-1 do begin
    if (Mag_B(i) lt B_cut) then L_B(i) = 10.0^((mb_sun - Mag_B(i)) / 2.5) $
    else L_B(i) = 0.0
    if (Mag_K(i) lt K_cut) then L_K(i) = 10.0^((mk_sun - Mag_K(i)) / 2.5) $
    else L_K(i) = 0.0
    if (Mag_B(i) lt B_cut and Mag_V(i) lt V_cut) then begin
        Mag_bj(i) = Mag_B(i) - 0.28 * (Mag_B(i) - Mag_V(i))
        L_bj(i) = 10.0^((mbj_sun - Mag_bj(i)) / 2.5)
    endif else begin
        Mag_bj(i) = 99.0
        L_bj(i) = 0.0
    endelse
endfor

;-------------------------------------------------------------------------------

; Create plot

; Open file for all plots
!p.charsize = 0.9
!p.charthick = 2
!p.thick = 2
!x.thick = 2
!y.thick = 2
!p.multi = 0
;print, 'Loading simple colour table for line graphics...'
red =   [0,1,1,0,0,0,1,1,0.5,0,0,0,0.5,0.5]
green = [0,1,0,1,0,1,0,1,0,0.5,0,0.5,0,0.5]
blue =  [0,1,0,0,1,1,1,0,0,0,0.5,0.5,0.5,0]
tvlct, 255*red,255*green,255*blue
; white background; black foreground
!p.background=1
!p.color=0

if (keyword_set(ps)) then begin
    set_plot, 'PS'
    device, filename = 'plot.ps', xsize = 12, ysize = 11, /color, $
      xoffset=1, $
      yoffset=5
endif
  
  xmin=12.0
  xmax=9.0
  ymin=0.000001
  ymax=0.1


  sel=where(Mag_B lt -17.0)
  mass1=alog10(StellarMass[sel])

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
	 /ylog,xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(Mh!U2!N/M!D!9ng!3!N)', $
	 ytitle = '!4U!3 (h!U3!NMpc!U-3!N(log!D10!N(M!U-1!N))'

  bin=0.1
  hist=(histogram(mass1,locations=c,min=9.0,max=12.0,binsize=bin))
  oplot,c+alog10(Hubble_h)+(bin/2.0),hist/(volume*0.1), color=70
  print, 'Stellar Mass Function'

  PRINT, FORMAT='(f12,",")', c-alog10(Hubble_h)+(bin/2.0)
  PRINT, FORMAT='(e15,",")', hist/(volume*0.1) 
  ;oplot,c+alog10(Hubble_h)-0.125,hist/(volume*0.1), color=70,linestyle=1


;DLB07
  c=[8.89510, 8.99510, 9.09510, 9.19510, 9.29510, 9.39510, 9.49510, 9.59510, $
     9.69510, 9.79510, 9.89510, 9.99510, 10.0951, 10.1951, 10.2951, 10.3951, $
     10.4951, 10.5951, 10.6951, 10.7951, 10.8951, 10.9951, 11.0951, 11.1951, $
     11.2951, 11.3951, 11.4951, 11.5951, 11.6951, 11.7951]
  hist=[0.0672752, 0.0610017, 0.0553210, 0.0497897, 0.0450793, 0.0409899, 0.0376410, $
        0.0344244, 0.0290783, 0.0236126, 0.0203428, 0.0181535, 0.0170988, 0.0155181, $
        0.0133526, 0.0112546, 0.00939909, 0.00777953, 0.00609034, 0.00444129,        $
        0.00288727, 0.00172032, 0.000889651, 0.000424755, 0.000255590, 0.000119603,  $
        6.79936e-05, 2.82624e-05, 1.06496e-05, 4.50560e-06]
 
  oplot,c,hist, color=70,linestyle=2

;Baldry2008

  symbols, 1, 0.4
  Readcol,'./data/baldrynew.dat',obsmass,a,fi,err,errdown,errup

  ;oploterror, obsmass+alog10(Hubble_h^2),fi/(Hubble_h^3), $
               ;errup/(Hubble_h^3)-fi/(Hubble_h^3),color = 200, errcolor = 200,$ 
               ;psym = 8,HATLENGTH = 70.0,/hiba
  ;oploterror, obsmass+alog10(Hubble_h^2),fi/(Hubble_h^3), $
               ;errdown/(Hubble_h^3)-fi/(Hubble_h^3),color = 200, errcolor = 200, $
               ;psym = 8,HATLENGTH = 70.0,/loba



 ;close,3
 ;openw,3,'./data/a.dat'
 ;for i=0,32 do 
    ;printf, 3,obsmass[i]+alog10(Hubble_h), $
              ;((errup[i]+errdown[i])/(2.0*Hubble_h^3))*(volume*0.1), $
              ;((errup[i]-errdown[i])/2.0)/(Hubble_h^3)*(volume*0.1)
 ;for i=33,35 do 
    ;printf, 3,obsmass[i]+alog10(Hubble_h),fi[i]/(Hubble_h^3)*volume*0.1, $
              ;(err[i])/(Hubble_h^3)*(volume*0.1)
 ;close,3

  symbols, 1, 0.5
  Readcol,'./data/baldryfit.dat',obsmass,fi,err
 
  oploterror, obsmass+alog10(Hubble_h), fi/(244141*0.1), $
              (obsmass/obsmass)-1.0+0.05, err/(244141*0.1), $ 
              color = 200, errcolor = 200, psym = 8,HATLENGTH = 80.0

  ;loadct,2
  ;symbols, 22, 0.15,color=200 
  oploterror, obsmass+alog10(Hubble_h), fi/(244141*0.1),err/(244141*0.1), $
              color = 200, errcolor = 200, psym = 8,HATLENGTH = 80.0

  loadct,6
  ;oplot, obsmass+alog10(Hubble_h), fi/(244141*0.1)+err/(244141*0.1)+sqrt(fi)/(244141*0.1),psym=8
  ;oplot, obsmass+alog10(Hubble_h), fi/(244141*0.1)-err/(244141*0.1)-sqrt(fi)/(244141*0.1),psym=8
  ;oploterror, obsmass+alog10(Hubble_h), fi/(244141*0.1), $ 
               ;err/(244141*0.1)+sqrt(fi)/(244141*0.1), color = 200, errcolor = 200, $
               ;psym = 8,HATLENGTH = 80.0

  ;oploterror, obsmass+alog10(Hubble_h), fi/(244141*0.1), err/(244141*0.1), $
              ;color = 200, errcolor = 200, psym = 8,HATLENGTH = 80.0



;BALDRY 2008 FIT

   Mstar = 4.446*(1e10)
   alpha1 = -0.46+1.0
   alpha2 = -1.58+1.0
   phistar1 = 0.00426
   phistar2 = 0.00058
   MM =  (findgen(30000)/1000)
   M=(10^MM)
   LLF = ((phistar1*((M/Mstar)^(alpha1)))+(phistar2*((M/Mstar)^(alpha2))))* exp(-M/Mstar)
   ;oplot,Alog10(M)+alog10(Hubble_h^2),(LLF*Alog(10))/(Hubble_h^3),color=200
   ;oplot,Alog10(M),(LLF*alog(10)),color=200
   L1=(LLF)/(Hubble_h^3)



;Bell et al. 2003

  symbols, 30, 0.15
  Readcol,'./data/bell2003gbandfit.dat',massbell,fi,err
  oploterror, massbell-0.102, fi/(244141*0.1),err/(244141*0.1), $
              color = 120, errcolor = 120, psym = 8,HATLENGTH = 80.0



  symbols, 1, 0.5,color=200
  xyouts,11.75 , 0.04, 'Baldry et al. (2008)'
  oplot, [11.85,11.85], [0.047,0.047], psym=8, color = 200
  
  symbols, 30, 0.15,color=120
  xyouts,11.75 , 0.018, 'Bell et al. (2003)'
  oplot, [11.85,11.85], [0.021,0.021], psym=8, color = 120
  

  ;oplot, [10.8,11.0], [0.0000023,0.0000023], linestyle=1,color=70
  oplot, [10.8,11.0], [0.0000046,0.0000046], linestyle=8,color=70
  oplot, [10.8,10.875],  [0.00001,0.00001], linestyle=2,color=70
  oplot, [10.925,11.0],  [0.00001,0.00001], linestyle=2,color=70
  ;xyouts, 10.75, 0.000004, 'MCMC Best Fit'
  ;xyouts, 10.75, 0.000004, 'Satellite Disruption'
  xyouts, 10.75, 0.000004, 'Delucia z=0'

  ;xyouts, 10.75, 0.000002, 'MCMC with free x shift'
  xyouts, 10.75, 0.0000089, 'De Lucia & Blaizot (2007)'

if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
endif

;------------------------------------------------------------------------------


end
