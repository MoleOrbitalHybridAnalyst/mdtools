// ./pdb2xyz in.pdb out.xyz 
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

   pdb.writeXYZ(argv[2],"system");

   return 0;

}
