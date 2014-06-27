import numpy
import sqlite3
import read_lgal

folder = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm/treedata/"
firstfile = 0 # must be 0
lastfile = 0
lastsnap = 61
Mgadget2Msun = 1.e10
select_num = 20
sample_prefix = 100
samplefile_prefix ="cut_optimalsample_allz_nh"
# use sqlite in memory
global db
global cursor
db = sqlite3.connect(':memory:')
cursor = db.cursor()
# create a table

cursor.execute('''CREATE TABLE tree (
curhalonr INTEGER,
filenr INTEGER,
treenr INTEGER,
halonr INTEGER, 
snapnum INTEGER,
mass REAL)
''')

cursor.execute('''CREATE TABLE selected (
filenr INTEGER,
treenr INTEGER,
halonr INTEGER, 
snapnum INTEGER,
mass REAL)
''')

ntreesfile = numpy.zeros(lastfile+1-firstfile,dtype=numpy.int32)
for ifile in range(firstfile,lastfile+1):
    filename = folder+"/trees_%03d.%d"%(lastsnap,ifile)
    f = open(filename,"rb")
    ntreesfile[ifile] = numpy.fromfile(f,numpy.int32,1)[0]
    f.close()


global firsttreeinfile
firsttreeinfile = numpy.cumsum(ntreesfile) - ntreesfile


# read tree input
(nHalos,nTrees,ngalstree,output_trees) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)

# set up indices
firsthalointree = numpy.cumsum(ngalstree)-ngalstree
lasthalointree = numpy.cumsum(ngalstree)


# insert data to database
count_halo = 0
for j in range(nTrees):
    nh = ngalstree[j]
    for i in range(nh):
        if(output_trees[count_halo]['FirstHaloInFOFgroup'] == i):
            filenr = output_trees[count_halo]['FileNr']
            mass = numpy.log10(output_trees[count_halo]['M_Crit200']*Mgadget2Msun)
            snapnum = output_trees[count_halo]['SnapNum']
            cursor.execute("INSERT INTO tree(curhalonr,filenr,treenr,halonr,snapnum,mass) VALUES (?,?,?,?,?,?)",(count_halo,int(output_trees[count_halo]['FileNr']),j,i,int(snapnum),float(mass)))
        count_halo += 1


min_mass = 10.0
max_mass = 18.0
step_mass = 0.25
nSteps = int((max_mass-min_mass)/step_mass)
nhalosbin = numpy.zeros((nSteps,lastsnap+1),dtype=numpy.int32)
selectednhalosbin = numpy.zeros((nSteps,lastsnap+1),dtype=numpy.int32)
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
    filenr = this_tree[index]['FileNr']
    mass = numpy.log10(this_tree[index]['M_Crit200']*Mgadget2Msun)
    cursor.execute("INSERT INTO selected (filenr,treenr,halonr,snapnum,mass) VALUES (?,?,?,?,?)",(int(filenr),int(treenr),int(index),int(snap),float(mass)))
    db.commit()


for i in range(nSteps):
    low_m = min_mass + i*step_mass
    high_m = min_mass + (i+1)*step_mass
    cursor.execute("SELECT * FROM tree WHERE mass BETWEEN ? AND ? AND snapnum = ? ORDER BY RANDOM()",(low_m,high_m,lastsnap))
    all_rows = cursor.fetchall()

    for data in all_rows[0:20]:
        treenr = data[2]
        this_tree = output_trees[firsthalointree[treenr]:lasthalointree[treenr]]
        this_tree_index = data[0]-firsthalointree[treenr]
        treecrawler(this_tree_index,this_tree,treenr)
        

for i in reversed(range(lastsnap+1)):
    for j in range(nSteps):
        low_m = min_mass + j*step_mass
        high_m = min_mass + (j+1)*step_mass
        cursor.execute("SELECT * FROM selected WHERE mass BETWEEN ? AND ? AND snapnum = ?",(low_m,high_m,i))
        result_selected = cursor.fetchall()
        cursor.execute("SELECT * FROM tree WHERE mass BETWEEN ? AND ? AND snapnum = ?",(low_m,high_m,i))
        result_tree = cursor.fetchall()
        this_Nselect = len(result_selected)
        this_Nhalos = len(result_tree)
        print "Snap",i,"step",j,":",this_Nselect,this_Nhalos
        if this_Nselect < min(this_Nhalos,select_num):
            rq_num = min(this_Nhalos,select_num)-this_Nselect
            cond = []
            for data in result_selected:
                string = "(filenr != %d OR treenr != %d OR halonr != %d OR snapnum != %d)"%(data[0],data[1],data[2],data[3])
                cond.append(string)
            cond_str = " AND ".join(cond)
            query_str = "SELECT * FROM tree WHERE mass BETWEEN ? AND ? AND snapnum = ? AND %s ORDER BY RANDOM() LIMIT %d"%(cond_str,rq_num)
            print query_str
            cursor.execute(query_str,(low_m,high_m,i))
            result_add = cursor.fetchall()
            for add_data in result_add:
                treenr = add_data[2]
                this_tree = output_trees[firsthalointree[treenr]:lasthalointree[treenr]]
                this_tree_index = add_data[0]-firsthalointree[treenr]
                treecrawler(this_tree_index,this_tree,treenr)
            print "add some more halos"
            cursor.execute("SELECT * FROM selected WHERE mass BETWEEN ? AND ? AND snapnum = ?",(low_m,high_m,i))
            result_selected = cursor.fetchall()
            this_Nselect = len(result_selected)
            print "Snap",i,"step",j,":",this_Nselect,this_Nhalos
    # get total halos in snap (selected)
    cursor.execute("SELECT * FROM selected WHERE snapnum = %d"%(i))
    total_nhalos = len(cursor.fetchall())
    f = open("%s_%d%03d"%(samplefile_prefix,sample_prefix,i),"w")
    f.write("%d\n"%total_nhalos)
    for j in range(nSteps):
        low_m = min_mass + j*step_mass
        high_m = min_mass + (j+1)*step_mass
        cursor.execute("SELECT * FROM selected WHERE mass BETWEEN ? AND ? AND snapnum = ?",(low_m,high_m,i))
        result_selected = cursor.fetchall()
        cursor.execute("SELECT * FROM tree WHERE mass BETWEEN ? AND ? AND snapnum = ?",(low_m,high_m,i))
        result_tree = cursor.fetchall()
        this_Nselect = float(len(result_selected))
        this_Nhalos = float(len(result_tree))
        if(this_Nselect > 0):
            weight = this_Nhalos/this_Nselect
            for data in result_selected:
                u_id = data[0]*1000000000000L+data[1]*1000000L+data[2]
                f.write("%ld\t%d\t%d\t%f\n"%(int(u_id),int(data[1]),int(data[0]),float(weight)))
    f.close()
db.close()


