# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:45:43 2017

@author: pyh
"""

def Readvasprun():
    import numpy as np
    import re
    fid = open('vasprun.xml','rt')
    lc = np.mat(np.zeros((3,3)))
    fir = True
    while fir == True:
        line = fid.readline()
        if 'NSW' in line:
            line = re.findall('(?<=> ).*(?=<)',line)
            nstep = int(line[0])
        if 'POTIM' in line:
            line = re.findall('(?<=> ).*(?=<)',line)
            step = int(line[0])
        if '<atoms>' in line:
            line = re.findall('(?<=> ).*(?=<)',line)
            natom = int(line[0])            
        if 'basis' in line:
            for direct in range(3):
                line = fid.readline()
                line = line.replace('<v>','')
                line = line.replace('</v>','')
                lc[direct,:] = np.mat(line)
                fir = False                          
    posa=[None]*(nstep+1)
    pos=np.mat(np.zeros((natom,3)))
        
    step=0
    atom=0        
    while step <= nstep or not line:
        line = fid.readline()
        if 'name="positions"' in line:
            for atom in range(natom):
                line = fid.readline()
                line = line.replace('<v>','')
                line = line.replace('</v>','')
                pos[atom,:] = np.mat(line)
            posa[step]=pos
            step += 1
    fid.close()
    return step,nstep,natom,lc,posa

def readCONTCAR():
    import numpy as np
    lc=np.mat(np.zeros((3,3)));
    fid = open('CONTCAR','rt')
    line = fid.readline()
    line = fid.readline()
    for i in range(0,3):
        lc[i,:]=np.mat(fid.readline())
    line = fid.readline()
    natom=np.mat(fid.readline())
    number=natom.sum()
    ap=np.mat(np.zeros((number,3)))
    while not 'Direct' in line:
        line = str(fid.readline())
    for i in range(0,number):
        line = fid.readline()
        if 'T' in line:
            line=line.replace('T','');
        if 'F' in line:
            line=line.replace('F','')
        ap[i,:]=np.mat(line)
    fid.close()
    return lc,ap

def nearby(lc,arg1,arg2):
    import numpy as np
    new_pos=arg2
    dp1=arg1*lc
    dp2=arg2*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    old_dis=0
    while np.all(old_dis != dis):
        old_dis=dis
        for i in range(0,3):
            for j in range(0,3):
                for k in range(0,3):
                    dp2=(new_pos+np.mat([i-1,j-1,k-1]))*lc
                    n_dis=(dp1-dp2)*np.transpose(dp1-dp2)
                    if np.all(n_dis < dis):
                        new_pos=new_pos+np.mat([i-1,j-1,k-1])
                        dis=n_dis
        arg2=new_pos
    return arg2  

def bond_length(ap,lc,arg1,arg2):
    import numpy as np
    import vasp
    import math
    arg1=arg1-1
    arg2=arg2-1
    ap[arg2,:]=vasp.nearby(lc,ap[arg1,:],ap[arg2,:])
    dp1=ap[arg1,:]*lc
    dp2=ap[arg2,:]*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    dis=math.sqrt(dis)
    return dis   

def bond_angle(ap,lc,arg1,arg2,arg3):
    import numpy as np
    import vasp
    import math
    arg1=arg1-1
    arg2=arg2-1
    arg3=arg3-1
    ap[arg1,:]=vasp.nearby(lc,ap[arg2,:],ap[arg1,:])
    ap[arg3,:]=vasp.nearby(lc,ap[arg2,:],ap[arg3,:])
    dp1=ap[arg1,:]*lc
    dp2=ap[arg2,:]*lc
    dp3=ap[arg3,:]*lc    
    ll1=dp2-dp1
    ll2=dp2-dp3
    cita=math.acos(ll1*np.transpose(ll2)/math.sqrt((ll1*np.transpose(ll1))*(ll2*np.transpose(ll2))))*180/math.pi
    return cita

def velcal(lc,posa,step):
    import numpy as np
    import vasp
    vela = [0]*(len(posa)-2)
    for i in range(len(vela)):
        vela[i]=np.mat(np.zeros((len(posa[1]),3)))
    for i in range(len(posa)-2):
        for j in range(len(posa[1])):
            print(posa[i+2][j,:])
            posa[i+2][j,:]=vasp.nearby(lc,posa[i][j,:],posa[i+2][j,:])
            print(posa[i+2][j,:])
            vela[i][j,:]=(posa[i+2][j,:]-posa[i][j,:])*lc/(2*step)
    return vela

#def findmax
