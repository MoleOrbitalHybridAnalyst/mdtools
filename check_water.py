#!/usr/bin/env python2
from pypdb import *
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_pdb', help='input pdb to be checked')
    parser.add_argument('-l', help='distance criteria', default='1.4')
    return parser.parse_args()

if __name__=="__main__":
    args = parse()
    pdb = pdb_obj(args.input_pdb)
    dist2 = float(args.l) ** 2

    defo = pdb_def("chainid W and atomtype O")
    selo = pdb.select_atoms(defo)
    if len(selo) == 0: raise Exception("no water in the input pdb")
    for index in selo:
        if pdb.pbc_distance2(index,index+1) > dist2 or \
            pdb.pbc_distance2(index, index+2) > dist2:
            print "serial",index + 1

    defhyd = pdb_def("chainid W and atomtype O and resname H3O")
    selhyd = pdb.select_atoms(defhyd)
    if len(selhyd) == 0: raise Exception("no hydronium in the pdb")
    elif len(selhyd) != 1: raise Exception("more than 1 hydronium in the pdb")
    index = selhyd[0]
    if pdb.pbc_distance2(index,index+1) > dist2 or \
        pdb.pbc_distance2(index, index+2) > dist2 or \
        pdb.pbc_distance2(index, index+3) > dist2:
        print "serial",index + 1
