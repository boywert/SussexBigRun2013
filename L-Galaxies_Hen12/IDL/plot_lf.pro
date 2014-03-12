;=========================================================================
  
   PRO plot_lf,mag,magname=magname,$
               ps=ps,box=box,Hubble_h=Hubble_h,$
               min=min,max=max,binsize=binsize

;=========================================================================
;
;  Example procedure to create luminosity function from magnitudes
;
;-------------------------------------------------------------------------

; Inputs:
; mag=array of magnitudes

magname = "M!DUV=GFUV!n - 5log!D10!nh"
ps = ps


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

; h
  x = x - (5*ALOG10(Hubble_h))
  y = y/(Hubble_h^3)

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
  

; Plot
  plot, x, y, /ylog,$
        xstyle=1, xrange = [xmin,xmax], xtitle = magname,$
        ystyle=1, yrange = 10.^[ylogmin,ylogmax],$
        ytitle = '!4U!3 /(h!u3!nMpc!u-3!nmag!u-1!n)'

;if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
;endif

;------------------------------------------------------------------------------


end
