;=========================================================================
;
;  Script to read in L-galaxies snapshot data
;  and produce a scatter plot of sfr versus mass

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 12
;   FirstFile = 0
;   LastFile = 511
    FirstFile = 5
    LastFile = 5
;
;-------------------------------------------------------------------------

@'/export/virgo/MillGas/WMAP/data/500_sam/Guo10_uv/snap_template.pro'
;
; Define which files you want to read in
if SNAP eq 12 then DirName = '/export/virgo/MillGas/WMAP/data/500_sam/Guo10_uv/'
if SNAP eq 12 then FileName = 'SA_z8.55'
ModelName = DirName + FileName

;-------------------------------------------------------------------------------

; read in data

read_lgal,ModelName,Template,GalStruct $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

StellarMass = (GalStruct.DiskMass+GalStruct.BulgeMass) * 1.0e10
SFR = GalStruct.sfr

;-------------------------------------------------------------------------------

; Create plot

start_plot

; Open file for all plots
if PS then set_plot, 'PS'
if PS then device, filename = 'mag.ps', /color, $
                  xsize = 15, ysize = 10, $
                  xoffset=1,  yoffset=5
  
x=StellarMass
y=SFR

plot, x, y, /nodata, /xlog, psym=3, $
  xtitle = 'Stellar mass/h!u-1!nM!d!Mn!n', $
  xrange = [1e2,1e12],xstyle=1, $
  ytitle = 'SFR / M!d!mn!n yr!u-1!n', $
  yrange = [-25,-5],ystyle=1
oplot, x, y, color=2, psym=3

if PS then device, /close_file
if PS then set_plot,'x'

;------------------------------------------------------------------------------
