import os

os.mkdir('NVT')
os.symlink('../INCAR_NVT','./NVT/INCAR')
os.symlink('../KPOINTS','./NVT/KPOINTS')
os.symlink('../POTCAR','./NVT/POTCAR')
os.symlink('../POSCAR','./NVT/POSCAR')
if os.path.exists('vdw_kernel.bindat'):
    os.symlink('../vdw_kernel.bindat','./NVT/vdw_kernel.bindat')
