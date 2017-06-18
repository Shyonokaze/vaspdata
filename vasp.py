# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:45:43 2017

@author: pyh
"""


def readvasprun():
    import xml.etree.ElementTree as ET 
    import numpy as np
    tree = ET.parse('vasprun.xml')
    root = tree.getroot()
    
    basis = True
    n=0
    lc = np.mat(np.zeros((3,3)))
    child = root[0]
    
    while child in root and basis:
        child = root[n]
        n += 1
        for grandson in child:
            if grandson.tag == 'atoms':
                natom = int(grandson.text)
            for ggrandson in grandson:
                if 'NSW' in ggrandson.attrib.values():
                    NSW = int(ggrandson.text)
                if 'POTIM' in ggrandson.attrib.values():
                    POTIM = float(ggrandson.text)
                if 'basis' in ggrandson.attrib.values():
                    for l in range(3):
                        lc[l,:] = np.mat(ggrandson[l].text) 
                    basis = False
    atom=0
    step=0
    
    posa=[None]*NSW
    for i in range(NSW):
        posa[i]=np.mat(np.zeros((natom,3)))
    
    for child in root:
        for grandson in child:
            for ggrandson in grandson:
                if 'positions' in ggrandson.attrib.values():
                    for atom in range(natom):
                        posa[step][atom,:]=np.mat(ggrandson[atom].text)
                    step += 1
    return POTIM,NSW,natom,lc,posa

def readpos(file): #此函数需针对selective dynamiscs的情况增加功能
    import numpy as np
    lc=np.mat(np.zeros((3,3)));
    fid = open(file,'rt')
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

def writepos(file,pos,vel=None):
    fido = open(file,'rt')
    fidn = open('new_POSCAR','wt')
    while True:
        line=fido.readline()
        print(line[0:(len(line)-1)],file=fidn)
        if 'Direct' in line:
            break
    fido.close()
    for i in range(len(pos)):
        for j in range(3):
            print('%.16f' % pos[i,j],file=fidn,end=' ') 
        print('',file=fidn)
    print('',file=fidn)
    if vel != None:
        for i in range(len(pos)):
            for j in range(3):
                print('%.12e' % vel[i,j],file=fidn,end=' ') 
            print('',file=fidn)

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
    arg1 -=1
    arg2 -=1
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
    arg1 -=1
    arg2 -=1
    arg3 -=1
    ap[arg1,:]=vasp.nearby(lc,ap[arg2,:],ap[arg1,:])
    ap[arg3,:]=vasp.nearby(lc,ap[arg2,:],ap[arg3,:])
    dp1=ap[arg1,:]*lc
    dp2=ap[arg2,:]*lc
    dp3=ap[arg3,:]*lc    
    ll1=dp2-dp1
    ll2=dp2-dp3
    cita=math.acos(ll1*np.transpose(ll2)/math.sqrt((ll1*np.transpose(ll1))*(ll2*np.transpose(ll2))))*180/math.pi
    return cita

def velcal(lc,posa,potim): #须针对只对一个原子计算速度的情况增加功能
    import numpy as np
    import vasp
    vela = [0]*(len(posa)-2)
    for i in range(len(vela)):
        vela[i]=np.mat(np.zeros((len(posa[1]),3)))
    for i in range(len(posa)-2):
        print(i) #测试
        for j in range(len(posa[1])):
            posa[i+2][j,:]=vasp.nearby(lc,posa[i][j,:],posa[i+2][j,:])
            vela[i][j,:]=(posa[i+2][j,:]-posa[i][j,:])*lc/(2*potim)
    return vela

def findlimit(vela,num,direct):
    import numpy as np
    num -=1
    if direct == 'x' or direct == 'X':
        direct = 1
    elif direct == 'y' or direct == 'Y':
        direct = 2
    elif  direct == 'z' or direct == 'Z':
        direct = 3
    direct -=1
    fvmax=0
    fvmin=0
    vmax=float(vela[0][num,direct])
    vmin=float(vela[0][num,direct])
    for i in range(len(vela)):
        if np.all(vmax < float(vela[i][num,direct])):
            vmax = float(vela[i][num,direct])
            fvmax = i
        if np.all(vmin > float(vela[i][num,direct])):
            vmax = float(vela[i][num,direct])
            fvmin = i
    return fvmax,fvmin,vela[fvmax],vela[fvmin]

def checkaway(lc,posa,vela,num):
    import numpy as np
    import math
    import vasp
    num -= 1
    nstep = len(posa)
    for n in range(7):
        step=int((n+1)*nstep/8)
        posa[step][num,:]=vasp.nearby(lc,posa[0][num,:],posa[step][num,:])
        dp=(posa[step][num,:]-posa[0][num,:])*lc
        dis = math.sqrt(dp*np.transpose(dp))
        if dis >8:
            check = True
            break
        elif dis > 5:
            dir=dp*np.transpose(vela[step][num,:])
            if dir > 0:
                check = True
                break
        check = False
    return check

#还需两个函数：二分法函数/创建新的目录及目录下的输入文件并移动到该目录下的函数
