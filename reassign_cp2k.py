#python3 reassign_cp2k.py input.pdb input.restart output.restart

#modify cp2k's restart file according to an input restart so that the configuration would be good-looking in vmd
#by the way, change center scheme to center setup_only

#requirements for the input pdb:
#   H3O must be together; the order must be OH2,H1,H2,H3
#   H2O must be together and ordered
#   atom type (element name) must be specified
#   chainid must be specified; water must have a chainid "W"
#   resname must be specified; H3O must have a resname "H3O"

#requirements for the input restart:
#   VELOCITY section must appear after COORD
#   CENTER appears only once 

from sys import argv
import re
import numpy as np

#extra acceptors; all the indexes starting from zero
#serca system:
#sero=[11801,11865,13907,12359]  #oxygen index of ser or ser-like residues
#sero=[11801,13907,12359]  #oxygen index of ser or ser-like residues
sero=[13907,12359]  #oxygen index of ser or ser-like residues
#sero=[13907]  #oxygen index of ser or ser-like residues
hisn=[14565]  #nritogen index of his which indeed accepts proton
his_aux=[14570]  #nitrogen index of his which is used for distance searching but does not accept proton actually
#serh=[11802,11866,13908,12360]
#serh=[11802,13908,12360]
serh=[13908,12360]
#serh=[13908]
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

#def dij(i,j):   #calulate distance between atom i and j; nopbc considered
#    tmp=(xs[i]-xs[j])**2+(ys[i]-ys[j])**2+(zs[i]-zs[j])**2
#    return np.sqrt(tmp)

def dij(i,j):   #calulate distance between atom i and j; nopbc considered
    tmp=(rxs[i]-rxs[j])**2+(rys[i]-rys[j])**2+(rzs[i]-rzs[j])**2
    return np.sqrt(tmp)

def switch(i,j): #switch coordinates between atom i and j
    tmp=xs[i]; xs[i]=xs[j]; xs[j]=tmp
    tmp=ys[i]; ys[i]=ys[j]; ys[j]=tmp
    tmp=zs[i]; zs[i]=zs[j]; zs[j]=tmp
    tmp=rxs[i]; rxs[i]=rxs[j]; rxs[j]=tmp
    tmp=rys[i]; rys[i]=rys[j]; rys[j]=tmp
    tmp=rzs[i]; rzs[i]=rzs[j]; rzs[j]=tmp
    tmp=vxs[i]; vxs[i]=vxs[j]; vxs[j]=tmp
    tmp=vys[i]; vys[i]=vys[j]; vys[j]=tmp
    tmp=vzs[i]; vzs[i]=vzs[j]; vzs[j]=tmp
    print("switch %d(%s) and %d(%s)"%(i+1,atomnames[i],j+1,atomnames[j]))

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

    #read restart to get the line number or "CENTER" and get the "COORD" and "VELOCITY" sections ; all the numbers starts from 1
    center_num=-1; coord_num=-1; velocity_num=-1; issubsys=0; line_num=0
    ratomnames=[]; rxs=[]; rys=[]; rzs=[]; rsegnames=[]; rresids=[]
    vxs=[]; vys=[]; vzs=[]
    fp_input_restart=open(argv[3],"r")
    for line in fp_input_restart:
        line_num+=1
        if re.match("\s+CENTER",line)!=None:
            center_num=line_num
        if re.match("\s*&SUBSYS",line)!=None: issubsys=1
        if issubsys==1 and re.match("\s*&END\s+SUBSYS",line)!=None: issubsys=0
        if issubsys==1: #if within the SUBSYS section
            if re.match("\s*&COORD",line)!=None: coord_num=line_num
            if re.match("\s*&END\s+COORD",line)!=None: 
                coord_num=-coord_num
                coord_end=line_num
            if coord_num>0 and line_num>coord_num:
                tmp=line.split()
                ratomnames.append(tmp[0])
                rxs.append(float(tmp[1]))
                rys.append(float(tmp[2]))
                rzs.append(float(tmp[3]))
                rsegnames.append(tmp[4])
                rresids.append(tmp[5])
            if coord_num!=-1 and re.match("\s*&VELOCITY",line)!=None:
                velocity_num=line_num
            if re.match("\s*&END\s+VELOCITY",line)!=None: 
                velocity_num=-velocity_num
                velocity_end=line_num
            if velocity_num>0 and velocity_num<line_num:
                tmp=line.split()
                vxs.append(float(tmp[0]))
                vys.append(float(tmp[1]))
                vzs.append(float(tmp[2]))
    fp_input_restart.close()
   # print(coord_num,velocity_num)
   # print(coord_end,velocity_end)


    #reassign until the donar H is selected again
    donarH_prev=np.array([1000,1000,1000]);donarH=-1;hyd_oindex=-1
    ####
    ttcount=0
    ####
    #while(np.sum((donarH-donarH_prev)**2)>0.01):
    while(ttcount<9):
        ttcount += 1
        #finding H3O; the index of H3O oxygen does not change 
        if hyd_oindex==-1:
            for index,(resname,atomtype) in enumerate(zip(resnames,atomtypes)):
                if resname=="H3O" and atomtype=="O":
                    hyd_oindex=index
                    #print(hyd_oindex)
                    break
        #finding the longest O-H bond in H3O
        hyd_hindex=-1; longest_bond=0.0
        for i in range(hyd_oindex+1,hyd_oindex+4):   #loop for H3O hydrogens
            #print(i)
            dtmp=dij(hyd_oindex,i)
            if(dtmp>longest_bond):
                hyd_hindex=i
                longest_bond=dtmp
                print("largest H3O bond = %f"%longest_bond)
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
            print(smallest_dist)
            #if min(dij(hish[num],hyd_hindex),dij(his_aux[num],hyd_hindex))> smallest_dist:
            if dij(hish[num],hisn[num])>smallest_dist:
                switch(hish[num],hyd_hindex)
        else: #ser is the acceptor
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
        
    #output new restart file
    coord_num=-coord_num; velocity_num=-velocity_num
    fp_input_restart=open(argv[3],"r")
    fp_output_restart=open(argv[4],"w")
    countline=0; isqm_kind=0; isqm_link=0; isqmmm=0; isqm_cell=0; isconstr=0
    for line in fp_input_restart:
        countline+=1
        if re.match("\s*&QMMM",line)!=None: isqmmm=1
        if re.match("\s*&END\s+QMMM",line)!=None: 
            isqmmm=0
            print("\n     @include qm_kind_cell_link\n",file=fp_output_restart)
        if re.match("\s*&CONSTRAINT",line)!=None: isconstr=1
        if re.match("\s*&END\s+CONSTRAINT",line)!=None: 
            isconstr=0
            #print("\n     @include constraint.inp\n",file=fp_output_restart)
            continue
        if re.match("\s*&QM_KIND",line)!=None: isqm_kind=1
        if re.match("\s*&END\s+QM_KIND",line)!=None: isqm_kind=0; continue
        if re.match("\s*&LINK",line)!=None: isqm_link=1
        if re.match("\s*&END\s+LINK",line)!=None: isqm_link=0; continue
        if re.match("\s*&MOTION",line)!=None: 
            print(" &MOTION\n     @INCLUDE constraint.inp",file=fp_output_restart)
            continue
        if isconstr==1: continue
        if isqmmm==1 and re.match("\s*&CELL",line)!=None: isqm_cell=1
        if isqmmm==1 and re.match("\s*&END\s+CELL",line)!=None: isqm_cell=0; continue
        if isqm_kind==1 or isqm_link==1 or isqm_cell==1: continue
        if countline==center_num: 
            print("     CENTER SETUP_ONLY",file=fp_output_restart)
        elif  countline==coord_num:
            print("     &COORD",file=fp_output_restart)
            for (atomname,x,y,z,segname,resid) in zip(ratomnames,rxs,rys,rzs,rsegnames,rresids):
                print("%s%26.16E%26.16E%26.16E %s %s"%(atomname,x,y,z,segname,resid),file=fp_output_restart)
        elif  countline==velocity_num:
            print("     &VELOCITY",file=fp_output_restart)
            for (vx,vy,vz) in zip(vxs,vys,vzs):
                print("       %26.16E%26.16E%26.16E"%(vx,vy,vz),file=fp_output_restart)
        elif  countline>coord_num and countline<coord_end:
            continue
        elif  countline>velocity_num and countline<velocity_end:
            continue
        else: print(line,file=fp_output_restart,end="")
    fp_input_restart.close()
    fp_output_restart.close()

    #write new pdb
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
    fp_pdb_inp.close()
    fp_pdb_out.close()
