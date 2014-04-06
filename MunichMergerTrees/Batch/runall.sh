#!/bin/bash

SimDir="/mnt/lustre/scratch/cs390/SUSSING2013_DATA/raw_subfind/"
SnapDirBase="snapdir_"
SnapFileBase="62.5_dm"
FirstSnapNr=0
LastSnapNr=61
SnapSkipFac=1
NumOutputFiles=8
FirstOutputFile=0
LastOutputFile=7
SnapFormat=1
FilesPerSnapshot=16
FileWithSnapList="/mnt/lustre/scratch/cs390/SUSSING2013_DATA/raw_subfind/zlist_MR7.txt"
BoxSize=62.5
Hashbits=8
BaseTree_param="B-BaseTree.param"


echo "OutputDir"  $SimDir > $BaseTree_param
echo "SnapshotFileBase"  $SnapFileBase >> $BaseTree_param
echo "SnapSkipFac"  $SnapSkipFac  >> $BaseTree_param
echo "LastSnapShotNr"  $LastSnapNr >> $BaseTree_param
echo "SnapFormat" $SnapFormat >> $BaseTree_param


HaloTree_param="B-HaloTree.param"
echo "OutputDir"  $SimDir > $HaloTree_param
echo "SnapshotFileBase"  $SnapFileBase >> $HaloTree_param
echo "FirstSnapShotNr"  $FirstSnapNr >> $HaloTree_param
echo "LastSnapShotNr"  $LastSnapNr >> $HaloTree_param
echo "SnapSkipFac"  $SnapSkipFac  >> $HaloTree_param
echo "NumberOfOutputFiles"  $NumOutputFiles  >> $HaloTree_param

TreeAddIDTab_param="L-TreeAddIDTab.param"
echo "SimulationDir"  $SimDir > $TreeAddIDTab_param
echo "LastSnapShotNr"  $LastSnapNr >> $TreeAddIDTab_param
echo "FirstFile"  $FirstOutputFile >> $TreeAddIDTab_param
echo "LastFile"  $LastOutputFile >> $TreeAddIDTab_param
echo "FilesPerSnapshot"  $FilesPerSnapshot  >> $TreeAddIDTab_param

TreeAddPosTab_param="L-TreeAddPosTab.param"
echo "SimulationDir"  $SimDir > $TreeAddPosTab_param
echo "LastSnapShotNr"  $LastSnapNr >> $TreeAddPosTab_param
echo "FirstFile"   $FirstOutputFile  >> $TreeAddPosTab_param
echo "LastFile"   $LastOutputFile  >> $TreeAddPosTab_param
echo "FilesPerSnapshot"  $FilesPerSnapshot  >> $TreeAddPosTab_param
echo "SnapshotFileBase"  $SnapFileBase  >> $TreeAddPosTab_param
echo "SnapFormat"  $SnapFormat  >> $TreeAddPosTab_param

TreeMakeUniqueIDs_param="L-TreeMakeUniqueIDs.param"
echo "SimulationDir" $SimDir > $TreeMakeUniqueIDs_param
echo "LastSnapShotNr" $LastSnapNr >> $TreeMakeUniqueIDs_param
echo "FirstFile"  $FirstOutputFile >> $TreeMakeUniqueIDs_param
echo "LastFile"  $LastOutputFile >> $TreeMakeUniqueIDs_param
echo "FilesPerSnapshot"  $FilesPerSnapshot >> $TreeMakeUniqueIDs_param
echo "FileWithSnapList"  $FileWithSnapList >> $TreeMakeUniqueIDs_param
echo "BoxSize"  $BoxSize  >> $TreeMakeUniqueIDs_param
echo "Hashbits"  $Hashbits  >> $TreeMakeUniqueIDs_param


current_dir=$(pwd)

cd ../Codes/B-BaseTree/
make
cd $current_dir
for i in {0..61}
do
../Codes/B-BaseTree/B-BaseTree ./B-BaseTree.param $( printf '%03d' $i)
done


cd ../Codes/B-HaloTrees
make
cd $current_dir
../Codes/B-HaloTrees/B-HaloTrees ./B-HaloTree.param 

cd ../Codes/L-TreeAddIDTab/
make
cd $current_dir
../Codes/L-TreeAddIDTab/L-TreeAddIDTab ./L-TreeAddIDTab.param

cd ../Codes/L-TreeAddPosTab/ 
make
cd $current_dir
../Codes/L-TreeAddPosTab/L-TreeAddPosTab ./L-TreeAddPosTab.param

cd ../Codes/L-TreeMakeUniqueIDs/ 
make
cd $current_dir
../Codes/L-TreeMakeUniqueIDs/L-TreeMakeUniqueIDs ./L-TreeMakeUniqueIDs.param

 
