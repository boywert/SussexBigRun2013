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
cursor.execute('''
CREATE TABLE tree(index  INTEGER,
                      filenr INTEGER,  
                      treenr INTEGER, 
                      halonr INTEGER,  
                      mass REAL)
''')

count_halo = 0
for nh in ngalstree:
    for i in range(nh):
        if(output_trees[count_halo]['FirstHaloInFOFgroup'] == i):
            print output_trees[count_halo]
        count_halo += 1


db.close()
