#!/bin/env python2

from __future__ import print_function
from pypdb import *
from sys import argv,stderr

if __name__ == '__main__':

    if len(argv) < 2:
        print('no arguments', file = stderr); exit()
    elif len(argv) == 2:
        in_pdb = argv[1]; out_pdb = argv[1]
    elif len(argv) == 3:
        in_pdb = argv[1]; out_pdb = argv[2]
    else:
        print('too many arguments', file = stderr); exit()

    pdb = pdb_obj(in_pdb)
    water_indexes = pdb.select_atoms(pdb_def('resname SPCF'))
    for i in water_indexes: pdb.resnames[i] = 'TIP3'

    pdb.write2file(out_pdb)
