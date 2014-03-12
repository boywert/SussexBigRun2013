;=========================================================================
;  Script to read in L-galaxies snapshot data
@'/mnt/lustre/scratch/virgo/SAM_output/Hen12_uv/snap_template.pro'
;-------------------------------------------------------------------------
DirName =  '/mnt/lustre/scratch/virgo/SAM_output/Hen12_uv/snapdir_013/'
FileName = 'SA_z7.88'
ModelName = DirName + FileName
FirstFile = 1
LastFile = 511
;-------------------------------------------------------------------------------
; read in data
if n_elements(Galstruct) eq 0 then read_lgal,ModelName,Template,G $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian
;------------------------------------------------------------------------------
redshift='z7.88'
set_plot, 'PS'
device, filename = redshift+'.ps', xsize = 12, ysize = 11, /color, $
        xoffset=1, $
        yoffset=5

x1=(g.bulgemass+g.diskmass)*1.0e10 ;stellar mass
x2=(g.sfr) ;SFR (Msol/yr)
y1=(g.metalscoldgas/g.coldgas) ;metallicity from gas
y2=(g.metalsbulgemass+g.metalsdiskmass)/(g.bulgemass+g.diskmass) ;metallicity from stars

;plot, alog10(x1),y2,title=redshift , xtitle='log(M!D*!N/M!Dsol!N)', ytitle='Z!Dgas', psym=3, xrange=[8.5,11.5]
plot,alog10(x2),y2,title=redshift , xtitle='log(SFR(!DM!Dsol!N!N/yr!N)', ytitle='Z!Dstars', psym=3
device, /close_file
set_plot,'x'
