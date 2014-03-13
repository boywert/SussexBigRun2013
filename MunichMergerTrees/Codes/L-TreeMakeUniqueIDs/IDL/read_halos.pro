
Base = "/afs/mpa/sim/milliMill/"
BaseOut = "/afs/mpa/sim/users/volker/milliMill/"
Num = 63
Nfiles = 1 ;8


HaloStruct = {$        
               OldDescendant           : 0L, $
               OldFirstProgenitor      : 0L, $
               OldNextProgenitor       : 0L, $
               OldFirstHaloInFOFgroup  : 0L, $
               OldNextHaloInFOFgroup   : 0L, $
               Len                  : 0L, $
               M_Mean200            : 0.0, $
               M_Crit200            : 0.0, $ 
               M_TopHat             : 0.0, $
               Pos                  : fltarr(3), $
               Vel                  : fltarr(3), $
               VelDisp              : 0.0, $
               Vmax                 : 0.0, $
               Spin                 : fltarr(3), $
               MostBoundID          : lon64arr(1), $
               SnapNum              : 0L, $ 
               FileNr               : 0L, $
               SubhaloIndex         : 0L, $
               SubhalfMass          : 0.0 $
             }

HaloStructAux= {$        
               HaloID               : lon64arr(1), $
               FileTreeNr           : lon64arr(1), $
               FirstProgenitor      : lon64arr(1), $
               LastProgenitor       : lon64arr(1), $
               NextProgenitor       : lon64arr(1), $
               Descendant           : lon64arr(1), $
               FirstHaloInFOFgroup  : lon64arr(1), $
               NextHaloInFOFgroup   : lon64arr(1), $
               Redshift             : 0.0D,        $
               PeanoKey             : 0L,          $
               Dummy                : 0L           $
             }


for FileNr = 0, Nfiles-1 do begin
    print, "Loading Filenr=", FileNr

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname= base+"/treedata/trees_"+ exts + ". " + string(FileNr)
    fname=strcompress(fname,/remove_all)

    openr,1, fname
 
    fname= base+"/treedata/tree_dbids_"+ exts + ". " + string(FileNr)
    fname=strcompress(fname,/remove_all)

    openr,2, fname
 

    fname= baseout+"/halos_" + string(FileNr) +".txt"
    fname=strcompress(fname,/remove_all)

    openw,3, fname, width=10000
 


    Ntrees = 0L
    TotNHalos = 0L
    readu,1,Ntrees,TotNhalos
    TreeNHalos = lonarr(Ntrees)
    readu,1,TreeNhalos

    for tr=0L, Ntrees-1 do begin

        Halo = replicate(HaloStruct, TreeNhalos(tr))
        readu,1, Halo

        HaloAux = replicate(HaloStructAux, TreeNhalos(tr))
        readu,2, HaloAux


        for i=0L, TreeNhalos(tr)-1 do begin
            printf, 3, HaloAux(i).HaloID(0), HaloAux(i).FirstProgenitor(0), HaloAux(i).LastProgenitor(0), $
  HaloAux(i).NextProgenitor(0), Halo(i).Len, Halo(i).M_Mean200, Halo(i).Pos(0:2), HaloAux(i).Redshift, HaloAux(i).PeanoKey
        endfor

    endfor

    close,3
    close,2
    close,1

endfor


end





