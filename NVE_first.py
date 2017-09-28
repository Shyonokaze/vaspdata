import os
import shutil
import sys
sys.path.append('/public/home/users/ruc001/bin/aBest')
import vasp

num=int(sys.argv[1])
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

os.chdir('../')
