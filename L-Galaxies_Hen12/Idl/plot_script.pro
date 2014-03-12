;=========================================================================
;
;  Example script to produce a plot from L-Galaxies data
;  This assumes that you have already used read_script.pro
;  to read the data into an L-Galaxies structure named G

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        PS=0b will draw on teh screen
;        PS=1b will save instead to a postscript file
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 63
;    SNAP=41
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
if PS then set_plot, 'PS'
if PS then device, filename = 'gas.ps', /color, $
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

if PS then device, /close_file
if PS then set_plot,'x'

;------------------------------------------------------------------------------
