#!/bin/bash
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -m bea
#$ -j y
#$ -cwd
#$ -pe openmpi 8 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

cubep3m_boxsize="47"
cubep3m_mesh="3456"
cubep3m_node="6"

pid_flag="1"
buffer_size="1.5"
drho="200"
n_chunks_pd="3"
n_chunks_total="27"

workspace="/home/c/cs/cs390/SussexBigRun2013/AHF_halos/"
particle_folder="/research/prace/cubepm_131025_6_1728_47Mpc_ext2/results/"
ahf_exec="/home/c/cs/cs390/SussexBigRun2013/ahf-v1.0-056.SUSSEXBIGRUN/bin/AHF-v1.0-056"
chunk_folder="/mnt/lustre/scratch/cs390/tmp/cubepm_131025_6_1728_47Mpc_ext2/chunked_output/"
ahfoutput_folder="/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131025_6_1728_47Mpc_ext2/"


while read line
do
    redshift=$(printf '%3.3f' $line)
    #make folder prepared for chunking
    for i in $(seq 0 $n_chunks_total)
    do
	this_workspace=$(printf '%s/z_%s_%d/chunk_%d/' $workspace $redshift $drho $i)
	mkdir -p "$this_workspace"
	this_chunkfolder=$(printf '%s/z_%s/chunk_%d/' $chunk_folder $redshift $i)
	mkdir -p "$this_chunkspace"
    done
    this_chunk_param=$(printf '%s/z_%s_%d/chunk_param' $workspace $redshift $drho)
    echo ${redshift} > ${this_chunk_param}
    echo "dummy" >>  $this_chunk_param
    echo $particle_folder >> $this_chunk_param
    echo $chunk_folder >> $this_chunk_param
    echo $cubep3m_boxsize >> $this_chunk_param
    echo $cubep3m_node >> $this_chunk_param
    echo $cubep3m_mesh >> $this_chunk_param
    echo $pid_flag >> $this_chunk_param
    echo $buffer_size >> $this_chunk_param
    echo $n_chunks_pd >> $this_chunk_param
    echo $n_chunks_pd >> $this_chunk_param
    echo $n_chunks_pd >> $this_chunk_param
done < halofinds

##mpirun -np 8 ../bin/AHF-v1.0-056 AHF.input-template2
