#include <iostream>
#include <string>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "pdb/pdb.h"
#include "pdb/pdbdef.h"
using namespace std;
using namespace PDB_NS;

void print_help(ostream& s) {
   s << "Usage:\n";
   s << "\treorder -f input_pdb -o output_pdb -d atom.def\n";
}

void message_abort(ostream& s, string str) {
   s << str << "\n";
   print_help(s);
   abort();
}

int main(int argc,char **argv) {

   /*--- read arguments ---*/

   char *fvalue=NULL,*ovalue=NULL,*dvalue=NULL;
   char argu_key;

   opterr = 0;
   
   while((argu_key = getopt(argc, argv,"f:o:d:")) != -1)
      switch(argu_key) {
         case 'f':
            fvalue = optarg;
            break;
         case 'o':
            ovalue = optarg;
            break;
         case 'd':
            dvalue = optarg;
            break;
         case '?':
            if(optopt == 'f') 
               message_abort(cerr, "ERROR: Option -f requires an argument");
            else if(optopt == 'o') 
               message_abort(cerr, "ERROR: Option -o requires an argument");
            else if(optopt == 'd') 
               message_abort(cerr, "ERROR: Option -s requires an argument");
            else message_abort(cerr, "ERROR: Unknown option -"+ string(1,optopt));
            break;
         default:
            message_abort(cerr, "ERROR: Unknown error");
      }

   if(!fvalue) message_abort(cerr, "ERROR: Option -f not found");
   if(!ovalue) message_abort(cerr, "ERROR: Option -o not found");
   if(!dvalue) message_abort(cerr, "ERROR: Option -d not found");

   /*--- end of read arguments ---*/

   /*--- read pdb and pdbdef ---*/
   PDB pdb(fvalue);
   PDBDef def(dvalue);
   auto cen = pdb.geoCenter(def);
   auto sel = pdb.selectAtoms(def);
   cout << "atom selection: ";
   for(auto i : sel) {
      cout << i << ' ';
   }
   cout << endl << "original center: ";
   cout << cen[0] << ' ' << cen[1] << ' ' << cen[2] << endl;
   pdb.shiftToMiddle(def);
   //pdb.shiftToMiddle(sel);
   pdb.write2file(ovalue);
}
