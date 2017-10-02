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

NP=sys.argv[1]
PBS_NODEFILE=sys.argv[2]
num1=int(sys.argv[3])
num2=int(sys.argv[4])

run='/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np '+NP+' -machinefile '+PBS_NODEFILE+' --mca btl self,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout'
try:
    os.mkdir('NVT')
    os.chdir('NVT')
    os.symlink('../INCAR_NVT','./INCAR')
    os.symlink('../KPOINTS','./KPOINTS')
    os.symlink('../POTCAR','./POTCAR')
    os.symlink('../POSCAR','./POSCAR')

    os.system(run)

    potim,nsw,natom,lc,posa=vasp.readvasprun()

    vela=vasp.velcal(lc,posa,potim)
    velas=[[0,0,0] for i in range(len(vela))]

    for i in range(len(vela)):
        velas[i]=(vela[i][num1]+vela[i][num2])/2

    fvmax,fvmin=vasp.findlimit(velas)

    vasp.writepos('POSCAR','pos_vmax',posa[fvmax],vela[fvmax])
    vasp.writepos('POSCAR','pos_vmin',posa[fvmin],vela[fvmin])

    shutil.copy('pos_vmax','../pos_vmax')
    shutil.copy('pos_vmin','../pos_vmin')

    os.chdir('../')
except:
    pass

try:
    os.mkdir('NVEmin')
    os.chdir('NVEmin')
    os.symlink('../INCAR_NVE','./INCAR')
    os.symlink('../KPOINTS','./KPOINTS')
    os.symlink('../POTCAR','./POTCAR')
    shutil.copy('../pos_vmax','./pos_vmax')
    shutil.copy('../record','./record')
    os.chdir('..')
except:
    pass

os.chdir('NVEmin')
ap,avel=vasp.readpos('pos_vmax')
avel[num1-1,2]=avel[num1-1,2]-0.05
avel[num2-1,2]=avel[num2-1,2]-0.05
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
    vela1=vasp.velcal(lc,posa,potim,num1)
    vela2=vasp.velcal(lc,posa,potim,num2)
    check1=vasp.checkaway(lc,posa,vela,num1)
    check2=vasp.checkaway(lc,posa,vela,num2)

    ns=vasp.bisec(check1)

    if ns:
        ap,avel=vasp.readpos('POSCAR')
        avel[num1-1,2]=avel[num1-1,2]-ns/2
        avel[num2-1,2]=avel[num2-1,2]-ns/2
        vasp.writepos('POSCAR','pos_next',ap,avel)
    
    shutil.copy('./pos_next','../POSCAR')
    shutil.copy('./record','../record')
    
    fid = open('record','rt')
    line=fid.readlines()
    fid.close()
    
    os.chdir('../')
    
os.chdir('../')

try:
    os.mkdir('NVEmax')
    os.chdir('NVEmax')
    os.symlink('../INCAR_NVE','./INCAR')
    os.symlink('../KPOINTS','./KPOINTS')
    os.symlink('../POTCAR','./POTCAR')
    shutil.copy('../pos_vmin','./pos_vmin')
    shutil.copy('../record','./record')
    os.chdir('..')
except:
    pass

os.chdir('NVEmax')
ap,avel=vasp.readpos('pos_vmin')
avel[num1-1,2]=avel[num1-1,2]-0.05
avel[num2-1,2]=avel[num2-1,2]-0.05
vasp.writepos('pos_vmin','POSCAR',ap,avel)

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
    vela1=vasp.velcal(lc,posa,potim,num1)
    vela2=vasp.velcal(lc,posa,potim,num2)
    check1=vasp.checkaway(lc,posa,vela,num1)
    check2=vasp.checkaway(lc,posa,vela,num2)

    ns=vasp.bisec(check1)

    if ns:
        ap,avel=vasp.readpos('POSCAR')
        avel[num1-1,2]=avel[num1-1,2]-ns/2
        avel[num2-1,2]=avel[num2-1,2]-ns/2
        vasp.writepos('POSCAR','pos_next',ap,avel)
    
    shutil.copy('./pos_next','../POSCAR')
    shutil.copy('./record','../record')
    
    fid = open('record','rt')
    line=fid.readlines()
    fid.close()
    
    os.chdir('../')
