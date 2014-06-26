import numpy
import sqlite3
import read_lgal

folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/"
firstfile = 0
lastfile = 0
lastsnap = 61
Mgadget2Msun = 1.e10
(nHalos,nTrees,ngalstree,output_trees) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)
db = sqlite3.connect(':memory:')
cursor = db.cursor()
cursor.execute('''CREATE TABLE tree (
curhalonr INTEGER,
filenr INTEGER,
treenr INTEGER,
halonr INTEGER, 
mass REAL)
''')

count_halo = 0
for j in range(nTrees):
    nh = ngalstree[j]
    for i in range(nh):
        if(output_trees[count_halo]['SnapNum'] == lastsnap):
            filenr = output_trees[count_halo]['FileNr']
            mass = numpy.log10(output_trees[count_halo]['M_Crit200']*Mgadget2Msun)
            cursor.execute("INSERT INTO tree(curhalonr,filenr,treenr,halonr,mass) VALUES (?,?,?,?,?)",(count_halo,int(output_trees[count_halo]['FileNr']),j,i,float(mass)))
        count_halo += 1



cursor.execute('''SELECT * FROM tree ORDER BY RANDOM() LIMIT 20''')
all_rows = cursor.fetchall()

for data in all_rows:
    print data
db.close()
