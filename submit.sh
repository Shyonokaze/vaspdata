#!/bin/bash
#PBS -N pyh
#PBS -l nodes=2:ppn=24
#PBS -q snode

#OPENMPI
export OPEN_MPI=/public/software/mpi/openmpi/1.6.5/intel
export PATH=${OPEN_MPI}/bin:$PATH
export INCLUDE=${OPEN_MPI}/include:$INCLUDE
export LD_LIBRARY_PATH=${OPEN_MPI}/lib:${OPEN_MPI}/share:/public/software/compiler/intel/composer_xe_2015.2.164/compiler/lib/intel64:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`

num=88

python /public/home/users/ruc001/work/7-PYH/7-BlackP-MD/1-MD/2aBest/main.py $NP $PBS_NODEFILE $num
