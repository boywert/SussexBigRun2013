Num = 1023

Base= "/home/volker/Billennium/C02_400/"


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

sfrflag = 0
bhflag = 0

;--------------------------------------------------


skip = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/subhalo_ids_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    Offset = 0L

    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, Offset

    if fnr eq 0 then begin
        IDs = lonarr(TotNids)
    endif

    if Nids gt 0 then begin
        
        locIDs = lonarr(Nids)
        readu,1, locIDs

        IDs(skip:skip+Nids-1) = locIDs(*)
        skip+= Nids

    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "Total IDs in grps=", TotNids
print

end
