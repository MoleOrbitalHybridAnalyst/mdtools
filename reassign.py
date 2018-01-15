#modify a pdb so that the configuration would be good-looking in vmd
#requirements for the input:
#H3O must be together; the order must be OH2,H1,H2,H3
#H2O must be together and ordered
#atom type (element name) must be specified
#chainid must be specified; water must have a chainid "W"
#resname must be specified; H3O must have a resname "H3O"
from sys import argv
import re
import numpy as np

#extra acceptors; all the indexes starting from zero
#serca system:
sero=[11801,11865,13907,12359]  #oxygen index of ser or ser-like residues
hisn=[14565]  #nritogen index of his which indeed accepts proton
his_aux=[14570]  #nitrogen index of his which is used for distance searching but does not accept proton actually
serh=[11802,11866,13908,12360]
hish=[14566]

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

def dij(i,j):   #calulate distance between atom i and j; nopbc considered
    tmp=(xs[i]-xs[j])**2+(ys[i]-ys[j])**2+(zs[i]-zs[j])**2
    return np.sqrt(tmp)

def switch(i,j): #switch coordinates between atom i and j
    tmp=xs[i]; xs[i]=xs[j]; xs[j]=tmp
    tmp=ys[i]; ys[i]=ys[j]; ys[j]=tmp
    tmp=zs[i]; zs[i]=zs[j]; zs[j]=tmp

if __name__=="__main__":
    #read pdb, construct lists for holding the pdb
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")
    atomnames=[];resnames=[];chainids=[];resids=[];xs=[];ys=[];zs=[];occs=[];tempfs=[];segnames=[];atomtypes=[];
    for line in fp_pdb_inp:
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
    #end of reading pdb

    #reassign until the donar H is selected again
    donarH_prev=np.array([1000,1000,1000]);donarH=-1;hyd_oindex=-1
    while(np.sum((donarH-donarH_prev)**2)>0.01):
        #finding H3O; the index of H3O oxygen does not change 
        if hyd_oindex==-1:
            for index,(resname,atomtype) in enumerate(zip(resnames,atomtypes)):
                if resname=="H3O" and atomtype=="O":
                    hyd_oindex=index
                    break
        #finding the longest O-H bond in H3O
        hyd_hindex=-1; longest_bond=0.0
        for i in range(hyd_oindex+1,hyd_oindex+4):   #loop for H3O hydrogens
            dtmp=dij(hyd_oindex,i)
            if(dtmp>longest_bond):
                hyd_hindex=i
                longest_bond=dtmp
        if sum((donarH_prev-np.array([1000,1000,1000]))**2)<0.01:
            print("index of H3O oxygen:",hyd_oindex+1,"(%s)"%atomnames[hyd_oindex],"index of H3O hydrogen:",hyd_hindex+1,"(%s)"%atomnames[hyd_hindex])
        #finding the closest acceptor
        close_owindex=-1; smallest_dist=1000
        for index,(chainid,atomtype) in enumerate(zip(chainids,atomtypes)):
            if chainid=="W" and atomtype=="O":
                dtmp=dij(index,hyd_hindex)
                if dtmp<smallest_dist:
                    close_owindex=index; smallest_dist=dtmp
        close_serindex=-1; close_hisindex=-1
        for index in sero:
            dtmp=dij(index,hyd_hindex)
            if dtmp<smallest_dist:
                close_serindex=index; smallest_dist=dtmp
        for (nindex,auxindex) in zip(hisn,his_aux):
            dtmp=dij(nindex,hyd_hindex)
            if dtmp<smallest_dist:
                close_hisindex=nindex; smallest_dist=dtmp
            else:
                dtmp=dij(auxindex,hyd_hindex)
                if dtmp<smallest_dist:
                    close_hisindex=nindex; smallest_dist=dtmp
        if close_serindex==-1 and close_hisindex==-1: #water is the acc.
            close_aindex=close_owindex
            switch(hyd_oindex,close_owindex)
            counttmp=0
            for i in range(hyd_oindex+1,hyd_oindex+4):  
                if i==hyd_hindex: continue
                counttmp+=1
                switch(i,close_owindex+counttmp)
        elif close_hisindex!=-1: #his is the acceptor
            close_aindex=close_hisindex
            for num,index in enumerate(hisn):
                if close_hisindex==index: break
            if min(dij(hish[num],hyd_hindex),dij(his_aux[num],hyd_hindex))> smallest_dist:
                switch(hish[num],hyd_hindex)
        else:
            close_aindex=close_serindex
            for num,index in enumerate(sero):
                if close_serindex==index: break
            #No. num serine has the closest acceptor
            dtmp=dij(serh[num],hyd_hindex)
            if dtmp>smallest_dist:
                switch(serh[num],hyd_hindex)
        print("index of closest acceptor:",close_aindex+1,"(%s)"%atomnames[close_aindex],"distance:",smallest_dist)
        if sum((donarH_prev-np.array([1000,1000,1000]))**2)>0.01:
            donarH_prev=donarH
        else:
            donarH_prev=np.array([-1000,-1000,-1000])
        donarH=np.array([xs[hyd_hindex],ys[hyd_hindex],zs[hyd_hindex]])
        


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
                pass
                print("{:6s}***** {:^4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}".format("ATOM",atomnames[countline],resnames[countline],chainids[countline],resids[countline],xs[countline],ys[countline],zs[countline],occs[countline],tempfs[countline],segnames[countline],atomtypes[countline]),file=fp_pdb_out)
            else:
                pass
                print("{:6s}{:5d} {:^4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}".format("ATOM",countline+1, atomnames[countline],resnames[countline],chainids[countline],resids[countline],xs[countline],ys[countline],zs[countline],occs[countline],tempfs[countline],segnames[countline],atomtypes[countline]),file=fp_pdb_out)

