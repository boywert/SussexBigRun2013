set_plot, 'PS'
device, filename = 'steve.ps', xsize = 12, ysize = 11, /color, $
        xoffset=1, $
        yoffset=5

x=g.mag[5]
y=g.magdust[5]

a = y-x

yarray = FltArr(7)
;print, size(yarray)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-16), count1)]=0
var[where(g.mag[5] lt (-17), count2)]=0
;print, count1, count2, size(var,/N_ELEMENTS)
result = TOTAL(var)
yarray[0]= result/(size(var,/N_ELEMENTS) - count1 - count2)
;print, size(var,/N_ELEMENTS) - count1 - count2

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-17), count1)]=0
var[where(g.mag[5] lt (-18), count2)]=0
result = TOTAL(var)
yarray[1]= result/(size(var,/N_ELEMENTS) - count1 - count2)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-18), count1)]=0
var[where(g.mag[5] lt (-19), count2)]=0
result = TOTAL(var)
yarray[2]= result/(size(var,/N_ELEMENTS) - count1 - count2)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-19), count1)]=0
var[where(g.mag[5] lt (-20), count2)]=0
result = TOTAL(var)
yarray[3]= result/(size(var,/N_ELEMENTS) - count1 - count2)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-20), count1)]=0
var[where(g.mag[5] lt (-21), count2)]=0
result = TOTAL(var)
yarray[4]= result/(size(var,/N_ELEMENTS) - count1 - count2)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-21), count1)]=0
var[where(g.mag[5] lt (-22), count2)]=0
result = TOTAL(var)
yarray[5]= result/(size(var,/N_ELEMENTS) - count1 - count2)

var = g.magdust[5] - g.mag[5]
var[where(g.mag[5] gt (-22), count1)]=0
var[where(g.mag[5] lt (-23), count2)]=0
result = TOTAL(var)
yarray[6]= result/(size(var,/N_ELEMENTS) - count1 - count2)

xarray=[-16.5,-17.5,-18.5,-19.5,-20.5,-21.5,-22.5]

print, yarray



;plot, y,a, /psym,xrange=[-24,-16],  xtitle='x', ytitle='y'
plot, g.mag[5], g.magdust[5]-g.mag[5],xrange=[-24,-16],title='z=4.89' , xtitle='M!Dins', ytitle='Af!Duv', psym=3, transparency=50
oplot, xarray,yarray, psym=1, color=2
device, /close_file
set_plot,'x'