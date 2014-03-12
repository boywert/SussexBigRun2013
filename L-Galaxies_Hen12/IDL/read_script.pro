
;=========================================================================
;
;  Script to read in L-galaxies snapshot data

;@'/export/data/virgo/MR7/data/500_sam/Guo10_uv/snap_template.pro'

@'/mnt/lustre/scratch/virgo/SAM_output/Hen12_uv/snap_template.pro'
;
;-------------------------------------------------------------------------

; Define which files you want to read in

;DirName =  '/mnt/lustre/scratch/virgo/SAM_output/Hen12_uv/snapdir_019/'

DirName =  '/mnt/lustre/scratch/phys/sc558/output/MR/newmodelb/snapdir_019/'

FileName = 'SA_z4.89'
ModelName = DirName + FileName
FirstFile = 1
LastFile = 511

;-------------------------------------------------------------------------------

; read in data

if n_elements(Galstruct) eq 0 then read_lgal,ModelName,Template,G $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian
;-------------------------------------------------------------------------------