# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 12:18:36 2017

@author: pyh
"""

import sys
sys.path.append('/public/home/users/ruc001/bin')

import vasp

num=52

potim,nsw,natom,lc,posa=vasp.readvasprun()
vela=vasp.velcal(lc,posa,potim,num)
check=vasp.checkaway(lc,posa,vela,num)

ns=vasp.bisec(check)

if ns:
    ap,avel=vasp.readpos('POSCAR')
    avel[num,2]=avel[num,2]-next
    vasp.writepos('POSCAR','pos_next',ap,avel)
