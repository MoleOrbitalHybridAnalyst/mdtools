#! /home/chhli/miniconda3/bin/python3
#python check_water.py in.pdb
#requirements for the input:
    #atom type (element name) must be specified
    #H2O must be in the order OH2 H1 H2
    #H3O must be in the order OH3 H1 H2 H3
    #chainid must be specified; water must have a chainid "W"
    #box must be specified
from sys import argv
import re
import numpy as np

def ext(line,m1,m2):
    return re.match("^\s*(\S+)\s*",line[m1:m2]).group(1)

def fill4(string):
    #python fills spaces in the left when using {:^ns}, this func fixes that
    if len(string)==1:
        return " "+string+"  "
    elif len(string)==2:
        return " "+string+" "
    elif len(string)==3:
        return " "+string
    elif len(string)==4:
        return string

def delta(x,y,i):
    d=x-y
    if(d>box[i]/2):
        d-=box[i]
    else:
        if(d<-box[i]/2):
            d+=box[i]
    return d

def dij(i,j):   #calulate distance between atom i and j
    #tmp=(xs[i]-xs[j])**2+(ys[i]-ys[j])**2+(zs[i]-zs[j])**2
    tmp=delta(xs[i],xs[j],0)**2+delta(ys[i],ys[j],1)**2+delta(zs[i],zs[j],2)**2
    return np.sqrt(tmp)

if __name__=="__main__":
    #read pdb, construct lists for holding the pdb
    box=[0,0,0]
    if len(argv) == 3:
       rbond = float(argv[2])
    else:
       rbond = 1.35
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")
    atomnames=[];resnames=[];chainids=[];resids=[];xs=[];ys=[];zs=[];occs=[];tempfs=[];segnames=[];atomtypes=[];
    for line in fp_pdb_inp:
        if re.match("^CRYST",line)!=None:
            tmp = line.split()
            if len(tmp)<4:
                print("box not defined")
                exit()
            box[0] = float(tmp[1]) 
            box[1] = float(tmp[2]) 
            box[2] = float(tmp[3]) 
        #only read ATOM section
        if re.match("^ATOM",line)==None:
            continue
        else:
           atomnames.append(fill4(ext(line,12,16)))
           resnames.append(ext(line,17,21))
           chainids.append(ext(line,21,22))
           resids.append(int(ext(line,22,26)))
           xs.append(float(ext(line,30,38)))
           ys.append(float(ext(line,38,46)))
           zs.append(float(ext(line,46,54)))
           occs.append(float(ext(line,54,60)))
           tempfs.append(float(ext(line,60,66)))
           segnames.append(ext(line,72,76))
           atomtypes.append(ext(line,76,78))
    fp_pdb_inp.close()
    if len(box)!=3 or box[0]==0.0 or box[1]==0.0 or box[2]==0.0:
        print("box not assigned")
        exit()
    #end of reading pdb
    
    for index,(chainid,atomtype) in enumerate(zip(chainids,atomtypes)):
        if chainid=="W" and atomtype=="O":
            if dij(index,index+1) > rbond or dij(index,index+2)> rbond:
                print(index+1)
    for index,(resname,atomtype) in enumerate(zip(resnames,atomtypes)):
        if resname=="H3O" and atomtype=="O":
            if dij(index,index+1) > rbond or dij(index,index+2)>rbond \
                    or dij(index,index+3) > rbond:
                print(index+1)
