from pypdb import *
from sys import argv

pdb = pdb_obj(argv[1])

defo = pdb_def("atomtype O and chainid W")
defh = pdb_def("atomtype H and chainid W")
defhyd = pdb_def("resname H3O and atomtype O and chainid W")

pdb.reorder_water(defo, defh, defhyd)

if len(argv) > 2:
   pdb.write2file(argv[2])
