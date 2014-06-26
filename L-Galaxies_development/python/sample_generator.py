import numpy
import sqlite3
import read_lgal

folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/"
firstfile = 0
lastfile = 0
lastsnap = 61
Mgadget2Msun = 1.e10
select_num = 20
(nHalos,nTrees,ngalstree,output_trees) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)
firsthalointree = numpy.cumsum(ngalstree)-ngalstree
print firsthalointree
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

min_mass = 10.0
max_mass = 18.0
step_mass = 0.25
nSteps = int((max_mass-min_mass)/step_mass)

for i in range(nSteps):
    low_m = min_mass + i*step_mass
    high_m = min_mass + (i+1)*step_mass
    cursor.execute("SELECT * FROM tree WHERE mass BETWEEN ? AND ? ORDER BY RANDOM()",(low_m,high_m))
    all_rows = cursor.fetchall()
    
    if len(all_rows) >= select_num:
        print len(all_rows)
        for data in all_rows[0:20]:
            treenr = data[2]
            this_tree = output_trees[firsthalointree[treenr],ngalstree[treenr]]
            print this_tree
db.close()
