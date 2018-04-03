from pypdb import *
from sys import argv

pdb = pdb_obj(argv[1])
for i in range(pdb.natom):
    if pdb.atomtypes[i] == 'H':
        pdb.occs[i] = 1.0079
    elif pdb.atomtypes[i] == 'C':
        pdb.occs[i] = 12.0107
    elif pdb.atomtypes[i] == 'N':
        pdb.occs[i] = 14.0067
    elif pdb.atomtypes[i] == 'O':
        pdb.occs[i] = 15.9994
pdb.write2file(argv[1])

