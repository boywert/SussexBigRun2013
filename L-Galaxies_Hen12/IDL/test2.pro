renske6x =[-0.04,0.41,0.77,1.01,1.26,1.51,1.77,2.03]
renske6y =[0.01197,0.00426,0.00173,0.00110,0.00026,0.00014,0.00002,0.00002]
renske6yerr =[0.00262,0.00089,0.00037,0.00024,0.00008,0.00004,0.00002,0.00002]
plot, renske6x, renske6y, /ylog
oploterr, renske6x, renske6y, renske6yerr, psym=3
;plot, renske6x, renske6y


;openw,1,'testing.txt' ;open the file to write
;for i=1,6,1 do printf,1,renske6x[i],renske6y[i] ;write data to file
;close, 1 ;close file


