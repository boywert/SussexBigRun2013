#!/bin/bash
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -m bea
#$ -j y
#$ -cwd
#$ -pe openmpi 216 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

mpi_chunk=216
mpi_ahf=8

cubep3m_boxsize=47
cubep3m_mesh=3456
cubep3m_node=6

pid_flag=1
buffer_size=1.5
drho=200
n_chunks_pd=3
n_chunks_total=27

workspace="/home/c/cs/cs390/SussexBigRun2013/AHF_halos/"
particle_folder="/research/prace/cubepm_131025_6_1728_47Mpc_ext2/results/"
ahf_folder="/home/c/cs/cs390/SussexBigRun2013/ahf-v1.0-056.SUSSEXBIGRUN/"
ahf_exec="/home/c/cs/cs390/SussexBigRun2013/ahf-v1.0-056.SUSSEXBIGRUN/bin/AHF-v1.0-056"
chunk_srcfolder="/home/c/cs/cs390/SussexBigRun2013/Chunking/"
chunk_exec="/home/c/cs/cs390/SussexBigRun2013/Chunking/chunk"
chunk_folder="/mnt/lustre/scratch/cs390/tmp/cubepm_131025_6_1728_47Mpc_ext2/chunked_output/"
ahfoutput_folder="/mnt/lustre/scratch/cs390/AHF_halos/cubepm_131025_6_1728_47Mpc_ext2/"

#compile things
cd ${ahf_folder}
make clean
make 
cd ${chunk_srcfolder}
make clean
make 


while read line
do
    
    redshift=$(printf '%3.3f' $line)
    echo "redshift = " $redshift
    firstfile=$(printf '%s/%sxv0.dat' $particle_folder $redshift)
    #make folder prepared for chunking
    for i in $(seq 0 $n_chunks_total)
    do
	this_workspace=$(printf '%s/z_%s_%d/chunk_%d/' $workspace $redshift $drho $i)
	mkdir -p "$this_workspace"
	this_chunkfolder=$(printf '%s/z_%s/chunk_%d/' $chunk_folder $redshift $i)
	mkdir -p "$this_chunkfolder"
    done
    rm -rf ${this_chunkfolder}/*
    cd ${this_workspace}
    this_chunk_param="chunk_param"

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
    if [ -e $firstfile ] 
    then
	this_pbs="submit.pbs"
	chunk_job_name=$(printf 'chunking_%s' $redshift)
	echo "#!/bin/bash" > $this_pbs
	echo "#$ -N" $chunk_job_name >> $this_pbs
	echo "#$ -M cs390@sussex.ac.uk" >> $this_pbs
	echo "#$ -m bea" >> $this_pbs
	echo "#$ -j y" >> $this_pbs
	echo "#$ -cwd" >> $this_pbs
	echo "#$ -pe openmpi" $mpi_chunk >> $this_pbs 
	echo "#$ -q mps_amd.q" >> $this_pbs
	echo "#$ -S /bin/bash" >> $this_pbs
	echo "module add sge" >> $this_pbs
	echo 'mpirun -np' $mpi_chunk $chunk_exec $this_chunk_param >> $this_pbs
	cat $this_pbs
	qsub $this_pbs
    fi
done < halofinds
##mpirun -np 8 ../bin/AHF-v1.0-056 AHF.input-template2
