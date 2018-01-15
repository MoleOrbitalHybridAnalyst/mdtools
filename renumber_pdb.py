#assume chainname is specified in the pdb
#will renumber atom index in each chain
#also add element symbol according to atom name given
from sys import argv
import re

#dictionary for mapping from atom name to element symbol
dic={"SOD": "NA", "POT": "K", "CLA": "CL"}

if __name__=="__main__":
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")
    index=0
    chain_prev="idontthinkachainhasthisname"
    for line in fp_pdb_inp:
        #read chain name and assign index
        if re.match("^ATOM",line)==None:
            #if not an ATOM line, throw out directly
            print(line,end="")
        else:
            #if an ATOM line, read chainname first
            chain_cur=line[72:76]
            if chain_cur!=chain_prev:
                index=1
                chain_prev=chain_cur  
            else:
                index+=1
            line=line[:6]+"{:5d}".format(index)+line[11:]
            #add element name
            element=line[12:16]
            element=re.match("^\s*(.+)",element).group(1)
            if element in dic:
                line=line[:76]+"{:>2s}".format(dic[element])+line[78:]
            else:
                line=line[:76]+"{:>2s}".format(element[0])+line[78:]
            print(line,end="")
