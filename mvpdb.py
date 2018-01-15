#python3 input.pdb output.pdb mx my mz
#output=input+{mx my mz}
#pdbs are in angstrom; mx my mz are in bohr
from sys import argv
import re

c=0.529177 #1 bohr =c angstrom

if __name__=="__main__":
    fn_pdb_inp=argv[1]
    fp_pdb_inp=open(argv[1],"r")
    xs=[];ys=[];zs=[]
    for line in fp_pdb_inp:
        if re.match("^ATOM",line)==None:
            continue
        else:
            xs.append(float(line[30:38])+c*float(argv[3]))
            ys.append(float(line[38:46])+c*float(argv[4]))
            zs.append(float(line[46:54])+c*float(argv[5]))
    fp_pdb_inp.close()
    fp_pdb_inp=open(argv[1],"r")
    fp_pdb_out=open(argv[2],"w")
    index=0
    for line in fp_pdb_inp:
        if re.match("^ATOM",line)==None:
            print(line,file=fp_pdb_out,end="")
            continue
        print(line[:30]+"{:8.3f}{:8.3f}{:8.3f}".format(xs[index],ys[index],zs[index])+line[54:],end="",file=fp_pdb_out)
        index+=1
    fp_pdb_inp.close()
    fp_pdb_out.close()
