#!/bin/bash
#PBS -N pyh
#PBS -l nodes=2:ppn=24
#PBS -q snode

#OPENMPI
export OPEN_MPI=/public/software/mpi/openmpi/1.6.5/intel
export PATH=${OPEN_MPI}/bin:$PATH
export INCLUDE=${OPEN_MPI}/include:$INCLUDE
export LD_LIBRARY_PATH=${OPEN_MPI}/lib:${OPEN_MPI}/share:/public/software/compiler/intel/composer
_xe_2015.2.164/compiler/lib/intel64:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`

num=88

python /public/home/users/ruc001/bin/aBest/NVT_mk.py
cd NVT/
/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl 
self,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout

python /public/home/users/ruc001/bin/aBest/NVT.py $num
cp ./pos_vmax ../
cp ./pos_vmin ../
cd ..

python /public/home/users/ruc001/bin/aBest/NVE_first.py $num

cd NVEmin/
tlast=$(grep 'end' ./record)
while [[ $tlast == '' ]]
do
python /public/home/users/ruc001/bin/aBest/NVE_mk.py 
td=$(tail -n1 ./record)
cd $td
cp ../record ./
/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl 
self,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout
python /public/home/users/ruc001/bin/aBest/NVE.py $num
mv ./pos_next ../POSCAR
mv ./record ../
cd ..
tlast=$(grep 'end' ./record)
done

cd ..
python /public/home/users/ruc001/bin/aBest/NVE_second.py $num

cd NVEmax/
tlast=$(grep 'end' ./record)
while [[ $tlast == '' ]]
do
python /public/home/users/ruc001/bin/aBest/NVE_mk.py
td=$(tail -n1 ./record)
cd $td
cp ../record ./
/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl
self,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout
python /public/home/users/ruc001/bin/aBest/NVE.py $num
mv ./pos_next ../POSCAR
mv ./record ../
cd ..
tlast=$(grep 'end' ./record)
done
