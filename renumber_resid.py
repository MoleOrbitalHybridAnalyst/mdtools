#! /home/chhli/miniconda3/bin/python3
# renumber_resid.py in.pdb
#vmd writepdb does not change resid
#psfgen is buggy if resid is not assigned properly
from sys import argv
import re

if __name__=="__main__":
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")

#construct a list of resids here
    resid=[]
    for i in range(9999):
        resid.append(int(i/3)+1)
#end of constructing the resids

    i=0
    for line in fp_pdb_inp:
        #read chain name and assign index
        if re.match("^ATOM",line)==None:
            #if not an ATOM line, throw out directly
            print(line,end="")
        else:
            #if an ATOM line, read chainname first
            line=line[:22]+"{:4d}".format(resid[i])+line[26:]
            print(line,end="")
            i+=1
