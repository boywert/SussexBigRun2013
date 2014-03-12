; Script to plot Disk/bulge split
    PS=0b
;   Which files
    FirstFile = 511
    LastFile = 511
;   Volume of box sampled by files (only correct for full volume)
    Volume=500^3*(LastFile-FirstFile+1.)/512. ; Volume in (Mpc/h)^3
    FileName = 'SA_z0.12'

;-------------------------------------------------------------------------------

.comp read_lgal
start_plot
xrange=[3e6,3e10]

; Open file for PS plot
if PS then set_plot, 'PS' else window,0
if PS then device, filename = 'sfh.ps', /color, $
                   xsize = 15, ysize = 10, $
                   xoffset=1,  yoffset=5

;DirName = '~/millgas/data/500_sam/Guo10_scaled/snapdir_056/'
;@'~/millgas/data/500_sam/Guo10_scaled/snap_template.pro'
DirName = '~/millgas/data/500_sam/sam_sfh2/snapdir_056/'
@'~/millgas/data/500_sam/sam_sfh2/snap_template.pro'
ModelName = DirName + FileName
read_lgal,ModelName,Template,Gnew $
         ,FirstFile=FirstFile,LastFile=LastFile;,/swap_endian
Gnew=Gnew[where(Gnew.BulgeMass+Gnew.DiskMass gt 0.1)]
Ngals=n_elements(Gnew)
Nbin=Gnew[0].sfh_ibin+1
x=Gnew[0].sfh_time[0:Nbin-1]
dt=Gnew[0].sfh_dt[0:Nbin-1]
y=fltarr(Nbin)
for i = 0, Nbin-1 do y[i]=total(Gnew.sfh_BulgeMass[i]+Gnew.sfh_DiskMass[i])/dt[i]
y=y*1e10/Volume
xnew=fltarr(Nbin+1)
xnew[Nbin]=3e6
xnew[0:Nbin-1]=x+0.5*dt
create_hist,xnew,y,xx,yy
print,'xnew=',xnew
print,'ynew=',y
plot, x, y, /nodata, /xlog, /ylog, $
      xrange = [3e6,3e10], xstyle=1, $
      ;yrange = [ymin, ymax], ystyle=1, $
      xtitle = 'Lookback time/yr', $
      ytitle = 'SFR / M!d!Mn!n yr!u-1!n (h!u-1!n Mpc)!u3!n'
oplot, xx, yy, color=2
oplot, x, y, psym=2, color=2
;xyouts, 1e7, 0.025, 'May 12 SFH binning', color=2

;DirName = '~/millgas/data/500_sam/sam_sfh2/snapdir_056/'
;@'~/millgas/data/500_sam/sam_sfh2/snap_template.pro'
DirName = '~/millgas/data/500_sam/sam_sfh4/snapdir_056/'
@'~/millgas/data/500_sam/sam_sfh4/snap_template.pro'
ModelName = DirName + FileName
read_lgal,ModelName,Template,Gold $
         ,FirstFile=FirstFile,LastFile=LastFile;,/swap_endian
Gold=Gold[where(Gold.BulgeMass+Gold.DiskMass gt 0.1)]
Ngals=n_elements(Gold)
Nbin=Gold[0].sfh_ibin+1
x=Gold[0].sfh_time[0:Nbin-1]
dt=Gold[0].sfh_dt[0:Nbin-1]
y=fltarr(Nbin)
for i = 0, Nbin-1 do y[i]=total(Gold.sfh_BulgeMass[i]+Gold.sfh_DiskMass[i])/dt[i]
y=y*1e10/Volume
xnew=fltarr(Nbin+1)
xnew[Nbin]=3e6
xnew[0:Nbin-1]=x+0.5*dt
print,'xold=',xnew
print,'yold=',y
create_hist,xnew,y,xx,yy
help,xx
help,yy
oplot, xx, yy, color=4
oplot, x, y, psym=2, color=4
;xyouts, 1e7, 0.015, 'Sep 12 SFH binning', color=4

if PS then device, /close_file
if PS then set_plot,'x'

