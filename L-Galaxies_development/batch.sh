#!/bin/bash
#$ -S /bin/bash
umask 002

#$ -cwd
#$ -q mps.q
#$ -l vf=16G
#$ -pe openmpi 128
#$ -j y
#$ -N Lgalaxy_all
#$ -o logs/$JOB_NAME_$JOB_ID.log


inputfolder="inputs_batch/"
mkdir -p $inputfolder
template="input/input_47mpc_template.par"
exec=./L-Galaxies
maxmemsize=4000


# No reionization
OutputDir="/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/"
SimulationDir="/mnt/lustre/scratch/cs390/47Mpc/"
XfracDir="/mnt/lustre/scratch/cs390/47Mpc/RT/47Mpc_f2_0_306/results/"
ReionizationOn=0
z0=50.0
zr=50.0
mkdir -p $OutputDir
mkdir -p $inputfolder
filename="${inputfolder}/input_0"

echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
echo "OutputDir" $OutputDir >> $filename
echo "SimulationDir" $SimulationDir >> $filename
echo "ReionizationOn" $ReionizationOn >> $filename
echo "Reionization_z0" $z0 >> $filename
echo "Reionization_zr" $zr >> $filename
echo "XfracDir" $XfracDir >> $filename
echo "MaxMemSize" $maxmemsize >> $filename
cat $template >> $filename


make metadata
mkdir -p "${OutputDir}/inputs/"
cp $filename "${OutputDir}/inputs/"
cp python/LGalaxyStruct.py "${OutputDir}/inputs/"


mpirun -np $NSLOTS $exec $filename > "${filename}_log"


# Okamoto 2008

OutputDir="/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/"
SimulationDir="/mnt/lustre/scratch/cs390/47Mpc/"
XfracDir="/mnt/lustre/scratch/cs390/47Mpc/RT/47Mpc_f2_0_306/results/"
ReionizationOn=2
z0=8.0
zr=7.0
mkdir -p $OutputDir
mkdir -p $inputfolder
filename="${inputfolder}/input_2"


echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
echo "OutputDir" $OutputDir >> $filename
echo "SimulationDir" $SimulationDir >> $filename
echo "ReionizationOn" $ReionizationOn >> $filename
echo "Reionization_z0" $z0 >> $filename
echo "Reionization_zr" $zr >> $filename
echo "XfracDir" $XfracDir >> $filename
echo "MaxMemSize" $maxmemsize >> $filename
cat $template >> $filename

make metadata
mkdir -p "${OutputDir}/inputs/"
cp $filename "${OutputDir}/inputs/"
cp python/LGalaxyStruct.py "${OutputDir}/inputs/"


mpirun -np $NSLOTS $exec $filename > "${filename}_log"




# patchy I

OutputDir="/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_I/"
SimulationDir="/mnt/lustre/scratch/cs390/47Mpc/"
XfracDir="/mnt/lustre/scratch/cs390/47Mpc/RT/47Mpc_f2_0_306/results/"
ReionizationOn=3
z0=50.0
zr=50.0
mkdir -p $OutputDir
mkdir -p $inputfolder
filename="${inputfolder}/input_3"


echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
echo "OutputDir" $OutputDir >> $filename
echo "SimulationDir" $SimulationDir >> $filename
echo "ReionizationOn" $ReionizationOn >> $filename
echo "Reionization_z0" $z0 >> $filename
echo "Reionization_zr" $zr >> $filename
echo "XfracDir" $XfracDir >> $filename
echo "MaxMemSize" $maxmemsize >> $filename
cat $template >> $filename


make metadata
mkdir -p "${OutputDir}/inputs/"
cp $filename "${OutputDir}/inputs/"
cp python/LGalaxyStruct.py "${OutputDir}/inputs/"


mpirun -np $NSLOTS $exec $filename > "${filename}_log"


# patchy II

OutputDir="/mnt/lustre/scratch/cs390/47Mpc/outputs/patchy_reionization_II/"
SimulationDir="/mnt/lustre/scratch/cs390/47Mpc/"
XfracDir="/mnt/lustre/scratch/cs390/47Mpc/RT/47Mpc_f2_0_306/results/"
ReionizationOn=4
z0=50.0
zr=50.0
mkdir -p $OutputDir
mkdir -p $inputfolder
filename="${inputfolder}/input_4"


echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
echo "OutputDir" $OutputDir >> $filename
echo "SimulationDir" $SimulationDir >> $filename
echo "ReionizationOn" $ReionizationOn >> $filename
echo "Reionization_z0" $z0 >> $filename
echo "Reionization_zr" $zr >> $filename
echo "XfracDir" $XfracDir >> $filename
echo "MaxMemSize" $maxmemsize >> $filename
cat $template >> $filename


make metadata
mkdir -p "${OutputDir}/inputs/"
cp $filename "${OutputDir}/inputs/"
cp python/LGalaxyStruct.py "${OutputDir}/inputs/"

mpirun -np $NSLOTS $exec $filename > "${filename}_log"


