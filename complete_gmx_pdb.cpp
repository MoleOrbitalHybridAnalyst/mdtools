// ./complete_gmx_pdb in.pdb out.pdb
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

const int MAXNUMWAT = 9999;

int main(int argc, char **argv) 
{
   PDB pdb(argv[1]);
   PDBDef def_protein("/home/chhli/share/protein.def");
   PDBDef def_water("/home/chhli/share/water.def");
   PDBDef def_lipid("/home/chhli/share/lipid.def");
   PDBDef def_pot("/home/chhli/share/pot.def");
   PDBDef def_cla("/home/chhli/share/cla.def");
   PDBDef def_sod("/home/chhli/share/sod.def");
   string chainid_protein = "P", chainid_lipid = "O";
   string chainid_ion = "O", chainid_water = "W";

   //do the protein
   vector<size_t> indexes = pdb.selectAtoms(def_protein);
   pdb.setChainid(def_protein, chainid_protein);
   int pre_resid = -1; int number = 1;
   for(auto index : indexes) {
      if(pdb.getResid(index) < pre_resid) number++;
      pre_resid = pdb.getResid(index);
      pdb.setSegname(index, chainid_protein + to_string(number));
   }

   //do the lipid
   pdb.setChainid(def_lipid, chainid_lipid);
   pdb.setSegname(def_lipid, chainid_lipid + "1");

   //do the ions
   number = 2;
   int resid = 1;
   indexes = pdb.selectAtoms(def_pot);
   if(!indexes.empty()) {
      for(auto index : indexes) {
         pdb.setSegname(index, chainid_ion + to_string(number));
         pdb.setChainid(index, chainid_ion);
         pdb.setResid(index, resid++);
      }
      number++;
   }
   indexes = pdb.selectAtoms(def_sod); resid = 1;
   if(!indexes.empty()) {
      for(auto index : indexes) {
         pdb.setSegname(index, chainid_ion + to_string(number));
         pdb.setChainid(index, chainid_ion);
         pdb.setResid(index, resid++);
      }
      number++;
   }
   indexes = pdb.selectAtoms(def_cla); resid = 1;
   if(!indexes.empty()) {
      for(auto index : indexes) {
         pdb.setSegname(index, chainid_ion + to_string(number));
         pdb.setChainid(index, chainid_ion);
         pdb.setResid(index, resid++);
      }
      number++;
   }

   //do the water
   indexes = pdb.selectAtoms(def_water); 
   int count = 0; number = 1;
   for(auto index : indexes) {
      if(count == MAXNUMWAT) {
         count = 0; number++;
      }
      pdb.setChainid(index, chainid_water);
      pdb.setSegname(index, chainid_water + to_string(number));
      pdb.setResid(index, count / 3 + 1);
      count++;
   }

   //write output
   pdb.write2file(argv[2]);
}
