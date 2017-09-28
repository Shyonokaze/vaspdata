# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 10:55:45 2017

@author: pyh
"""

import shutil
import sys
sys.path.append('/public/home/users/ruc001/bin/aBest')
import vasp

num=int(sys.argv[1])
potim,nsw,natom,lc,posa=vasp.readvasprun()
vela=vasp.velcal(lc,posa,potim)
fvmax,fvmin=vasp.findlimit(vela,num)

vasp.writepos('POSCAR','pos_vmax',posa[fvmax],vela[fvmax])
vasp.writepos('POSCAR','pos_vmin',posa[fvmin],vela[fvmin])

shutil.copy('pos_vmax','../pos_vmax')
shutil.copy('pos_vmin','../pos_vmin')
