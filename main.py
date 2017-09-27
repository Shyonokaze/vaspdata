# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 14:58:01 2017

@author: pyh
"""
import os
import shutil
import sys
sys.path.append('/public/home/users/ruc001/bin')
import vasp

run='/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE --mca btlself,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout'
num=52

os.mkdir('NVT')
os.chdir('NVT')
os.symlink('../INCAR_NVT','./INCAR')
os.symlink('../KPOINTS','./KPOINTS')
os.symlink('../POTCAR','./POTCAR')
os.symlink('../POSCAR','./POSCAR')

os.system(run)

potim,nsw,natom,lc,posa=vasp.readvasprun()

vela=vasp.velcal(lc,posa,potim)
fvmax,fvmin=vasp.findlimit(vela,num=52)

vasp.writepos('POSCAR','pos_vmax',posa[fvmax],vela[fvmax])
vasp.writepos('POSCAR','pos_vmin',posa[fvmin],vela[fvmin])

shutil.copy('pos_vmax','../pos_vmax')
shutil.copy('pos_vmin','../pos_vmin')

os.chdir('../')

os.mkdir('NVEmin')
os.chdir('NVEmin')
os.symlink('../INCAR_NVE','./INCAR')
os.symlink('../KPOINTS','./KPOINTS')
os.symlink('../POTCAR','./POTCAR')
shutil.copy('../pos_vmax','./pos_vmax')
shutil.copy('../record','./record')

ap,avel=vasp.readpos('pos_vmax')
avel[num-1,2]=avel[num-1,2]-0.1
vasp.writepos('pos_vmax','POSCAR',ap,avel)

fid = open('record','rt')
line=fid.readlines()
fid.close()

while 'end' not in line:
    vasp.filemv('INCAR')
    ind=line[-1].replace('\n','')
    os.chdir(ind)
    shutil.copy('../record','./record')
    os.system(run)
    
    potim,nsw,natom,lc,posa=vasp.readvasprun()
    vela=vasp.velcal(lc,posa,potim,num)
    check=vasp.checkaway(lc,posa,vela,num)

    ns=vasp.bisec(check)

    if ns:
        ap,avel=vasp.readpos('POSCAR')
        avel[num-1,2]=avel[num-1,2]-ns
        vasp.writepos('POSCAR','pos_next',ap,avel)
    
    shutil.copy('./pos_next','../POSCAR')
    shutil.copy('./record','../record')
    
    fid = open('record','rt')
    line=fid.readlines()
    fid.close()
    
    os.chdir('../')
