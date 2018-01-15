#include <iostream>
#include <string>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "/home/chhli/share/pdb/pdb.h"
using namespace std;
using namespace PDB_NS;


int main(int argc,char **argv) {
   PDB pdb(argv[1]);
   pdb.guessAllAtomtypes();
   PDBDef def_spcf;
   def_spcf.pushBack(PDBField::resname,"SPCF");
   vector<size_t> water_indexes = pdb.selectAtoms(def_spcf);
   for(auto iter = water_indexes.begin(); iter != water_indexes.end(); ++iter){
      pdb.setResname(*iter, "TIP3");
   }
   pdb.write2file(argv[2]);
}
