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
   s << "\tpdbpbc -f input_pdb -o output_pdb -b lx,hx,ly,hy,lz,hz\n";
}

void message_abort(ostream& s, string str) {
   s << str << "\n";
   print_help(s);
   abort();
}

int main(int argc,char **argv) {

   /*--- read arguments ---*/

   char *fvalue=NULL,*ovalue=NULL,*bvalue=NULL;
   char argu_key;

   opterr = 0;
   
   while((argu_key = getopt(argc, argv,"f:o:b:")) != -1)
      switch(argu_key) {
         case 'f':
            fvalue = optarg;
            break;
         case 'o':
            ovalue = optarg;
            break;
         case 'b':
            bvalue = optarg;
            break;
         case '?':
            if(optopt == 'f') 
               message_abort(cerr, "ERROR: Option -f requires an argument");
            else if(optopt == 'o') 
               message_abort(cerr, "ERROR: Option -o requires an argument");
            else if(optopt == 'b') 
               message_abort(cerr, "ERROR: Option -b requires an argument");
            else message_abort(cerr, "ERROR: Unknown option -"+ string(1,optopt));
            break;
         default:
            message_abort(cerr, "ERROR: Unknown error");
      }

   if(!fvalue) message_abort(cerr, "ERROR: Option -f not found");
   if(!ovalue) message_abort(cerr, "ERROR: Option -o not found");
   if(!bvalue) message_abort(cerr, "ERROR: Option -b not found");

   /*--- end of read arguments ---*/

   /*--- read pdb and pdbdef ---*/
   PDB pdb(fvalue);
   float lx,hx, ly, hy, lz, hz;
   sscanf(bvalue, "%f,%f,%f,%f,%f,%f", &lx, &hx, &ly, &hy, &lz, &hz);
   pdb.pbcWrap(lx, hx, ly, hy, lz, hz);
   pdb.write2file(ovalue);

   return 0;
}
