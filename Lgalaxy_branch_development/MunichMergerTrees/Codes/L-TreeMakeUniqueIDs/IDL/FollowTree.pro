
Base = "/afs/mpa/sim/milliMill/"

HaloStruct = {$        
               Descendant           : 0L, $
               FirstProgenitor      : 0L, $
               NextProgenitor       : 0L, $
               FirstHaloInFOFgroup  : 0L, $
               NextHaloInFOFgroup   : 0L, $
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

HaloAux    = {$        
               HaloID               : lon64arr(1), $
               FileTreeNr           : lon64arr(1), $
               FirstProgenitor      : lon64arr(1), $
               LastProgenitor       : lon64arr(1), $
               NextProgenitor       : lon64arr(1), $
               Descendant           : lon64arr(1), $
               FirstHaloInFOFgroup  : lon64arr(1), $
               NextHaloInFOFgroup   : lon64arr(1), $
               Redshift             : 0.0D $
             }

Num = 63
Nfiles =8

seed=42L

AllHalos= 0L

for FileNr = 0, Nfiles-1 do begin
    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname= base+"/treedata/trees_"+ exts + ". " + string(FileNr)
    fname=strcompress(fname,/remove_all)

    openr,1, fname
    Ntrees = 0L
    TotNHalos = 0L
    readu,1,Ntrees,TotNhalos
    close,1

    AllHalos+= TotNhalos
endfor

TreeAll = replicate(HaloStruct, AllHalos)
TreeAux = replicate(HaloAux, AllHalos)


AllHalos = 0L

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
 
    Ntrees = 0L
    TotNHalos = 0L
    readu,1,Ntrees,TotNhalos
    TreeNHalos = lonarr(Ntrees)
    readu,1,TreeNhalos

    for tr=0L, Ntrees-1 do begin

        Tree = replicate(HaloStruct, TreeNhalos(tr))
        readu,1, Tree

        TreeAll(AllHalos:AllHalos+TreeNhalos(tr)-1) = Tree(*)

        Tree = replicate(HaloAux, TreeNhalos(tr))
        readu,2, Tree

        TreeAux(AllHalos:AllHalos+TreeNhalos(tr)-1) = Tree(*)

        AllHalos+= TreeNhalos(tr)
     endfor
    close,1
    close,2
endfor




AllHalos = 0L

for FileNr = 0, Nfiles-1 do begin

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname= base+"/treedata/trees_"+ exts + ". " + string(FileNr)
    fname=strcompress(fname,/remove_all)

    openr,1, fname
    Ntrees = 0L
    TotNHalos = 0L
    readu,1,Ntrees,TotNhalos
    TreeNHalos = lonarr(Ntrees)
    readu,1,TreeNhalos

    for tr=0L, Ntrees-1 do begin

        Tree = replicate(HaloStruct, TreeNhalos(tr))
        readu,1, Tree
 
        halonr = fix(randomu(seed)*TreeNhalos(tr))

                    
        stack = lonarr(2000)
        stackptr = -1
    
        stack(++stackptr) = halonr

        Msum = -double(Tree(halonr).Len)
        count = -1L

        while stackptr ge 0 do begin

           p = stack(stackptr--)

           Msum += double(Tree(p).Len)

;           if (halonr eq 31) and (p ne halonr) then begin
;           print, Tree(p).Len, TreeAux(AllHalos+p).HaloID - TreeAux(AllHalos+halonr).HaloID
;           endif


           count++

           q = Tree(p).FirstProgenitor
                        
           while q ge 0 do begin
               stack(++stackptr) = q
               q = Tree(q).NextProgenitor
           endwhile
                         
        endwhile

           ind = where((TreeAux.HaloID ge TreeAux(AllHalos+halonr).FirstProgenitor) and $
                       (TreeAux.HaloID le TreeAux(AllHalos+halonr).LastProgenitor))

        if ind(0) ne -1 then begin
        print, total(TreeAll(ind).Len, /double) - Msum, count, TreeAux(AllHalos+halonr).LastProgenitor - $
        TreeAux(AllHalos+halonr).FirstProgenitor +1        , halonr
        endif else begin
             print , "failed",  Tree(halonr).FirstProgenitor
        endelse     


       AllHalos+= TreeNhalos(tr)

    endfor
    close,1
endfor



ende:
        close,1
        close,2
end

