#python3 shift_wrapping.py input_wrapped.pdb lx ly lz shiftx shifty shiftz
from sys import argv
import re
import numpy as np

if __name__=="__main__":
    xs=[];ys=[];zs=[]
    xl=0;xh=0;yl=0;yh=0;zl=0;zh=0
    fp=open(argv[1],"r")
    count=-1
    for line in fp:
        if re.match("^ATOM",line)==None:
            continue
        count+=1
        x=float(line[30:38]);y=float(line[38:46]);z=float(line[46:54])
        xs.append(x); ys.append(y); zs.append(z)
        if count==0:
            xl=x;xh=x;yl=y;yh=y;zl=z;zh=z
        if x>xh:
            xh=x
        if x<xl:
            xl=x
        if y>yh:
            yh=y
        if y<yl:
            yl=y
        if z>zh:
            zh=z
        if z<zl:
            zl=z
    fp.close()
    centerx=(xh+xl)/2+float(argv[5])
    centery=(yh+yl)/2+float(argv[6])
    centerz=(zh+zl)/2+float(argv[7])
    lx=float(argv[2])
    ly=float(argv[3])
    lz=float(argv[4])
    fp=open(argv[1],"r")
    count=-1
    for line in fp:
        if re.match("^ATOM",line)==None:
            print(line,end="")
            continue
        count+=1
        x,y,z=[xs[count],ys[count],zs[count]]
        if x<centerx-lx/2:
            x+=lx
        elif x>centerx+lx/2:
            x-=lx
        if y<centery-ly/2:
            y+=ly
        elif y>centery+ly/2:
            y-=ly
        if z<centerz-lz/2:
            z+=lz
        elif z>centerz+lz/2:
            z-=lz
        print(line[:30]+"{:8.3f}{:8.3f}{:8.3f}".format(x,y,z)+line[54:],end="")
    fp.close() 
