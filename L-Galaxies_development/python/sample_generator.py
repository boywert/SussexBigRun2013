import numpy
import sqlite3
import read_lgal

folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/"
firstfile = 0 # must be 0
lastfile = 0
lastsnap = 61
Mgadget2Msun = 1.e10
select_num = 20

# set up indices
firsthalointree = numpy.cumsum(ngalstree)-ngalstree
lasthalointree = numpy.cumsum(ngalstree)

# use sqlite in memory
db = sqlite3.connect(':memory:')
cursor = db.cursor()
# create a table

cursor.execute('''CREATE TABLE tree (
curhalonr INTEGER,
filenr INTEGER,
treenr INTEGER,
halonr INTEGER, 
mass REAL)
''')

ntreesfile = []
for ifile in range(firstfile,lastfile+1):
    filename = folder+"/trees_%03d.%d"%(lastsnap,ifile)
    f = open(filename,"rb")
    ntreesfile.append(numpy.fromfile(f,numpy.int32,1)[0])
    f.close()

# read tree input
(nHalos,nTrees,ngalstree,output_trees) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)

global firsttreeinfile
firsttreeinfile = numpy.cumsum(ntreesfile) - ntreesfile

# insert data to database
count_halo = 0
for j in range(nTrees):
    nh = ngalstree[j]
    for i in range(nh):
        if(output_trees[count_halo]['SnapNum'] == lastsnap) & (output_trees[count_halo]['FirstHaloInFOFgroup'] == i):
            filenr = output_trees[count_halo]['FileNr']
            mass = numpy.log10(output_trees[count_halo]['M_Crit200']*Mgadget2Msun)
            cursor.execute("INSERT INTO tree(curhalonr,filenr,treenr,halonr,mass) VALUES (?,?,?,?,?)",(count_halo,int(output_trees[count_halo]['FileNr']),j,i,float(mass)))
        count_halo += 1


min_mass = 10.0
max_mass = 18.0
step_mass = 0.25
nSteps = int((max_mass-min_mass)/step_mass)

# open files

global listsample
global listindex
global listtree
listsample = []
listindex = []
listtree = []

for i in range(lastsnap+1):
    listsample.append([])
    listindex.append([])
    listtree.append([])

def treecrawler(index,this_tree,treenr):
    if this_tree[index]['NextProgenitor'] > -1:
        treecrawler(this_tree[index]['NextProgenitor'],this_tree,treenr)
    if this_tree[index]['FirstProgenitor'] > -1:
        treecrawler(this_tree[index]['FirstProgenitor'],this_tree,treenr)
    snap = this_tree[index]['SnapNum']
    listsample[snap].append(this_tree[index])
    listindex[snap].append(index)
    listtree[snap].append(treenr-firsttreeinfile[this_tree[index]['SnapNum']])


for i in range(nSteps):
    low_m = min_mass + i*step_mass
    high_m = min_mass + (i+1)*step_mass
    cursor.execute("SELECT * FROM tree WHERE mass BETWEEN ? AND ? ORDER BY RANDOM()",(low_m,high_m))
    all_rows = cursor.fetchall()
    
    if len(all_rows) >= select_num:
        for data in all_rows[0:20]:
            treenr = data[2]
            this_tree = output_trees[firsthalointree[treenr]:lasthalointree[treenr]]
            this_tree_index = data[0]-firsthalointree[treenr]
            treecrawler(this_tree_index,this_tree)
db.close()

f = []
for i in range(lastsnap+1):
    f.append(open("sample.%03d"%(i),"w"))

for i in range(lastsnap+1):
    print len(listsample[i])
    for j in range(len(listsample[i])):
        print listindex[i][j],listtree[i][j]
for i in range(lastsnap+1):
    f[i].close()
