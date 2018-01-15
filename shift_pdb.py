#python3 shift.py input.pdb  shiftx shifty shiftz
#do pure shifting
from sys import argv
import re
import numpy as np

if __name__=="__main__":
    xs=[];ys=[];zs=[]
    fp=open(argv[1],"r")
    count=-1
    for line in fp:
        if re.match("^ATOM",line)==None:
            continue
        count+=1
        x=float(line[30:38]);y=float(line[38:46]);z=float(line[46:54])
        xs.append(x); ys.append(y); zs.append(z)
    fp.close()
    shiftx=float(argv[2])
    shifty=float(argv[3])
    shiftz=float(argv[4])
    fp=open(argv[1],"r")
    count=-1
    for line in fp:
        if re.match("^ATOM",line)==None:
            print(line,end="")
            continue
        count+=1
        x,y,z=[xs[count],ys[count],zs[count]]
        x+=shiftx
        y+=shifty
        z+=shiftz
        print(line[:30]+"{:8.3f}{:8.3f}{:8.3f}".format(x,y,z)+line[54:],end="")
    fp.close() 
