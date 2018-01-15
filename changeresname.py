#change resname according to the dict
from sys import argv
import re

#dictionary for mapping from the old resname to the new resname 
dic={"TIP3": "SPCF"}

if __name__=="__main__":
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")
    for line in fp_pdb_inp:
        if re.match("^ATOM",line)==None:
            print(line,end="")
        else:
            resname=line[17:21]
            resname=re.match("(.+)\s*",resname).group(0)
            if resname in dic:
                line=line[:17]+"{:4s}".format(dic[resname])+line[21:]
            print(line,end="")
    fp_pdb_inp.close()
