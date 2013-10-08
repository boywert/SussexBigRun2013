pro plot_procedure,Galstruct,ps=ps,snap=snap

;  Example procedure to produce a plot from L-Galaxies data.
;  This assumes that you have already used read_script.pro
;  to read the data into an L-Galaxies structure.
;  A producure is useful when you want to used begin..end blocks

;  Inputs:
;   Galstruct - L-Galaxies structure
;   ps - (optional) if set then write to file; otherwise to screen
;   snap - (optional) used to set plot range
;
  if n_elements(snap) eq 0 then snap=61
;
;-------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array

; Total stellar mass in Msun/h
  StellarMass = (GalStruct.DiskMass+GalStruct.BulgeMass) * 1.0e10
  ColdGas = GalStruct.ColdGas * 1.0e10
  HotGas = GalStruct.HotGas * 1.0e10

;-------------------------------------------------------------------------------

; Create plot

  start_plot

; Open file for all plots
  if keyword_set(ps) then set_plot, 'PS'
  if keyword_set(ps) then device, filename = 'gas.ps', /color, $
                                  xsize = 15, ysize = 10, $
                                  xoffset=1,  yoffset=5
  
  x=ColdGas
  y=Hotgas

  plot, x, y, /nodata, /xlog, /ylog, psym=3, $
        xtitle = 'Cold Gas mass/h!u-1!nM!d!Mn!n', $
        xrange = [1e6,1e12],xstyle=1, $
        ytitle = 'Hot Gas mass/h!u-1!nM!d!Mn!n', $
        yrange = [1e2,1e14],ystyle=1
  oplot, x, y, color=2, psym=3

  if keyword_set(ps) then device, /close_file
  if keyword_set(ps) then set_plot,'x'

end
