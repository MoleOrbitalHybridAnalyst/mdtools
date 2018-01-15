#python reorder.py in.pdb out.pdb
#modify a pdb so that the configuration would be good-looking in vmd
#requirements for the input:
    #atom type (element name) must be specified
    #H2O must be in the order OH2 H1 H2
    #H3O must be in the order OH3 H1 H2 H3
    #chainid must be specified; water must have a chainid "W"
    #box must be specified
#when extra_h is provided then a residue protonated config will be gen'ed
#else a config containing a hydronium will be gen'ed
from sys import argv
import re
import numpy as np

#extra_h = [19]   #indexes of extra H; only one extra_h is supported yet
extra_h = []
extra_o = []

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

def switch1d(x,y):
    return y,x

def switch(i,j): #switch coordinates between atom i and j
    tmp=xs[i]; xs[i]=xs[j]; xs[j]=tmp
    tmp=ys[i]; ys[i]=ys[j]; ys[j]=tmp
    tmp=zs[i]; zs[i]=zs[j]; zs[j]=tmp

if __name__=="__main__":
    #read pdb, construct lists for holding the pdb
    box=[0,0,0]
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
    print("reading finished")
    #end of reading pdb
    
    oxygens=[]      #indexes of the oxygens
    hydrogens=[]    #indexes of the hydrogens
    for index,(chainid,atomtype) in enumerate(zip(chainids,atomtypes)):
        if chainid=="W" and atomtype=="O":
            oxygens.append(index)
        if chainid=="W" and atomtype=="H":
            hydrogens.append(index)
    oxygens.extend(extra_o)
    hydrogens.extend(extra_h)
    destination = -1
    if len(extra_h)==0:     # if there are no extra_h, H3O should exist
        for index,(resname,atomtype) in enumerate(zip(resnames,atomtypes)):
            if resname=="H3O" and atomtype=="O":
                destination = index
                break
        destination += 3
    else:
        destination = extra_h[0]
    if destination == -1:
        print("no H3O found")
        exit()
    if 2*len(oxygens)+1!=len(hydrogens):
        print("2*len(oxygens)+1!=len(hydrogens)")
        exit()
    hydro_back = hydrogens[:]
    for io in oxygens:
        countj = 0; countj1 = -1; countj2 = -1
        for countj,jh in enumerate(hydrogens): #loop for H to find the bonded
            if jh==-1:
                continue
            if countj1==-1:
                d1 = dij(io,jh)
                countj1 = countj
            elif countj2==-1:
                d2 = dij(io,jh)
                countj2 = countj
                if d1 > d2:
                    d1,d2=switch1d(d1,d2)
                    countj1,countj2=switch1d(countj1,countj2)
            else:
                d2_ = dij(io,jh)
                if d2_ < d2 and d2_ > d1:
                    d2 = d2_
                    countj2 = countj
                elif d2_ < d1 and d2_ > d2:
                    d1 = d2_
                    countj1 = countj
                elif d2_ < d1 and d2_ < d2:
                    if d1 < d2:
                        d2 = d2_
                        countj2 = countj
                    else:
                        d1 = d2_
                        countj1 = countj
        # exclude the H just assigned
        flag = 0
        for countj,jh in enumerate(hydro_back):
            if jh==io+1:
                flag += 1
                hydrogens[countj] = -1
            if jh==io+2:
                flag += 1
                hydrogens[countj] = -1
        if flag!=2:
            print("cannot find two hydrogens of serial %d"%(io+1))
            exit()
        switch(io+1,hydro_back[countj1])
        switch(io+2,hydro_back[countj2])
    #all the water should be fine now
    #@@@
    #print(xs[799],ys[799],zs[799])
    #print(xs[800],ys[800],zs[800])
    #print(xs[801],ys[801],zs[801])
    #@@@
    flag = 0; excess_index = -1
    for jh in hydrogens:
        if jh!=-1:
            flag += 1
            excess_index = jh
    if flag!=1:
        print("number of excess protons != 1, why?")
        exit()
    if destination!=excess_index:
        print("the excess proton is not the third H of H3O or extra_o, why?")
        exit()
    d1 = box[0]+box[1]+box[2]; final_oxygen = -1
    for io in oxygens:
        d1_ = dij(io,excess_index)
        if d1_ < d1:
            d1 = d1_
            final_oxygen = io
    #@@@
    #print("excess_index = %d"%(excess_index+1))
    #print("destination = %d"%(destination+1))
    #@@@
    switch(destination,excess_index)
    if len(extra_h)==0: # if there is H3O in pdb
        switch(final_oxygen, destination - 3)
        switch(final_oxygen + 1, destination - 2)
        switch(final_oxygen + 2, destination - 1)
    else:
        print("%d"%(final_oxygen+1))

    # write pdb
    countline=-1
    fp_pdb_inp=open(argv[1],"r")
    fp_pdb_out=open(argv[2],"w")
    for line in fp_pdb_inp:
        #only read ATOM section
        if re.match("^ATOM",line)==None:
            print(line,end="",file=fp_pdb_out)
        else: 
            countline+=1
            if countline >= 99999:
                print("{:6s}***** {:^4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}".format("ATOM",atomnames[countline],resnames[countline],chainids[countline],resids[countline],xs[countline],ys[countline],zs[countline],occs[countline],tempfs[countline],segnames[countline],atomtypes[countline]),file=fp_pdb_out)
            else:
                print("{:6s}{:5d} {:^4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}".format("ATOM",countline+1, atomnames[countline],resnames[countline],chainids[countline],resids[countline],xs[countline],ys[countline],zs[countline],occs[countline],tempfs[countline],segnames[countline],atomtypes[countline]),file=fp_pdb_out)

