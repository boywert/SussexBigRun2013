import numpy
import sqlite3
import read_lgal

folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/"
firstfile = 0
lastfile = 0
lastsnap = 61
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
        if(output_trees[count_halo]['FirstHaloInFOFgroup'] == i):
            cursor.execute("INSERT INTO tree(curhalonr,filenr,treenr,halonr,mass) VALUES (?,?,?,?,?)",(count_halo,output_trees[count_halo]["FileNr"],j,i,2.))

#(count_halo,output_trees[count_halo]['FileNr'],j,i,output_trees[count_halo]['M_Crit200']))
            print output_trees[count_halo]
        count_halo += 1


db.close()
