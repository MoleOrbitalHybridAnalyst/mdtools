// ./pdb2ndx in.pdb out.ndx 
// ~/share/*def needed to determine groups
// assume no H3O as this is gmx
#include <iostream>
#include <string>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "pdb.h"
using namespace std;
using namespace PDB_NS;

int main(int argc, char **argv) 
{
   PDB pdb(argv[1]);
   PDBDef def_protein("/home/chhli/share/protein.def");
   PDBDef def_water("/home/chhli/share/water.def");
   PDBDef def_lipid("/home/chhli/share/lipid.def");
   //PDBDef def_pot("/home/chhli/share/pot.def");
   //PDBDef def_cla("/home/chhli/share/cla.def");
   //PDBDef def_sod("/home/chhli/share/sod.def");

   vector<Group> grps;
   grps.emplace_back("PROTEIN", pdb.selectAtoms(def_protein));
   grps.emplace_back("WATER", pdb.selectAtoms(def_water));
   grps.emplace_back("LIPID", pdb.selectAtoms(def_lipid));

   pdb.writeIndexFile(argv[2],grps);

   return 0;

}
