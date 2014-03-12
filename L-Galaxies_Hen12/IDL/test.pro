x=intarr[5]
print, x
x=[1,2,3,4,5,1,2,3]
print, x
y=[2,4,6,8,10,12,14,16]
z=[0,0,0,0,0]
a = [10,100,1000,10000,100000]
b = [5,50,500,5000,50000,500000]


for i=0,4 do begin
   y = x^2
   ;if (x[i] LT 3) then begin
   ;   x[i]=0


z = where(y LT 10, count)

print, x, y
print, count
print, "space"
print, y[z]
print, x[z]
