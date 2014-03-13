
Num = 1021

Base= "/gpfs/mpa/vrs/Billennium/halo_C02/C02_400/"

SnapSkipFac = 2

;--------------------------------------------------

if Num ge 1000 then begin
   exts='0000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
endif else begin
   exts='000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
endelse



f = Base + "/treedata/sub_desc_sf" + strcompress(string(SnapSkipFac),/remove_all)+"_"+exts

openr,1,f

TotNsubhalos = 0L

readu,1,TotNsubhalos
print, "TotNsubhalos=", TotNsubhalos, " in snap=",exts


descendant_halonr = lonarr(TotNsubhalos)
descendant_snapnr = lonarr(TotNsubhalos)

readu,1, descendant_halonr
readu,1, descendant_snapnr

close,1



end
