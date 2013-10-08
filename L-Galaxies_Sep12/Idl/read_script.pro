;=========================================================================
;
;  Script to read in L-galaxies snapshot data

; Read in the Template file describing the galaxy structure
@'/export/data/virgo/MillGas/WMAP7/data/500_sam/Guo10/snap_template.pro'
;@'~/millgas/data/500_sam/Guo10_sfh4/snap_template.pro'
;@'~/millgas/data/500_sam/Hen12_sfh4/snap_template.pro'
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName ='/export/data/virgo/MillGas/WMAP7/data/500_sam/Guo10/snapdir_061/'
;DirName = '~/millgas/data/500_sam/Hen12_sfh4/snapdir_061/'
;DirName = '../output/'
FileName = 'SA_z0.00'
ModelName = DirName + FileName
FirstFile = 511
LastFile = 511

;-------------------------------------------------------------------------------

; read in data

if n_elements(Galstruct) eq 0 then read_lgal,ModelName,Template,Galstruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------
