

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


Num = 63

 max_m200 = 0.0

for FileNr = 0, 511 do begin

print, "Doing FileNr= ", FileNr

Base = "/ptmp/vrs/Millennium/"


exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)


fname= base+"/treedata/trees_"+ exts + ". " + string(FileNr)
fname=strcompress(fname,/remove_all)

openr,1, fname
Ntrees = 0L
TotNHalos = 0L
readu,1,Ntrees,TotNhalos
print,"Ntrees= ", Ntrees
print,"TotNhalos= ", TotNhalos
TreeNHalos = lonarr(Ntrees)
readu,1,TreeNhalos

for tr=0L, Ntrees-1 do begin

    Tree = replicate(HaloStruct, TreeNhalos(tr))

    readu,1, Tree

    ind = where(Tree.SnapNum eq 18)

    if ind(0) ne -1 then begin	

       max_m200_local = max(Tree(ind).M_Mean200)

       if max_m200_local gt max_m200 then begin

           ind2 = where(max_m200_local eq Tree(ind).M_Mean200)

	   k = ind(ind2(0))
	   
	   Max_FileNr = FileNr
	   Max_TreeNr = Tr
           Max_Halo = k

	   print, "FileNr= ", FileNr,"  TreeNr= ", Tr, "  MaxHalo=", k, "   M200= ",  Tree(k).M_Mean200
           
           max_m200 = Tree(k).M_Mean200
       endif
    endif

endfor  
close,1

endfor

end
