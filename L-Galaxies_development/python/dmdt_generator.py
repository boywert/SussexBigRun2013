import numpy
import read_lgal

folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
firstfile = 0
lastfile = 127
lastsnap = 75

(nHalos,nTrees,ngalstree,treedata) = read_lgal.read_lgalinput(folder,firstfile,lastfile,lastsnap)

print tree
