pro create_hist,x,y,xx,yy,diff=diff

; Given points x defining edges of bins, and values y in bins, will
; return points to create a histogram
; Optionally will differentiate y across the interval

nx=n_elements(x)
ny=n_elements(y)
if nx ne ny+1 then message,'n_elements(x) != n_elements(y)+1'

nbin=2*nx
xx=fltarr(nbin)
yy=fltarr(nbin)

for ibin=0,nx-1 do begin
   xx[2*ibin]=x[ibin]
   xx[2*ibin+1]=x[ibin]
endfor

yy[0]=1e-30
for ibin=0,nx-2 do begin
   yy[2*ibin+1]=y[ibin]
   yy[2*ibin+2]=y[ibin]
endfor
yy[nbin-1]=1e-30

end

