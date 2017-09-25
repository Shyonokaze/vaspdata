# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 10:55:45 2017

@author: pyh
"""
import sys
sys.path.append('/public/home/users/ruc001/bin')

import vasp

potim,nsw,natom,lc,posa=vasp.readvasprun()

vela=vasp.velcal(lc,posa,potim)
fvmax,fvmin=vasp.findlimit(vela,num=52)

vasp.writepos('POSCAR','pos_vmax',posa[fvmax],vela[fvmax])
vasp.writepos('POSCAR','pos_vmin',posa[fvmax],vela[fvmax])
