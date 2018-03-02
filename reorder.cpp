#include <iostream>
#include <string>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "pdb/pdb.h"
using namespace std;
using namespace PDB_NS;

void print_help(ostream& s) {
   s << "Usage:\n";
   s << "\treorder -f input_pdb -o output_pdb -s segid4water [-g] [-a]\n";
}

void message_abort(ostream& s, string str) {
   s << str << "\n";
   print_help(s);
   abort();
}

//FILE* open_rfile(char *fname) {
//   if(!fname) {
//      print_help(cerr); abort();
//   } else {
//      FILE* fp = fopen(fname,"r");
//      if(!fp) {
//         string fname_str(fname);
//         message_abort(cerr, "ERROR: Cannot find file "+fname_str);
//      }
//      return fp;
//   }
//}
//
//
//FILE* open_wfile(char *fname) {
//   if(!fname) {
//      print_help(cerr); abort();
//   } else {
//      FILE* fp = fopen(fname,"w");
//      if(!fp) {
//         string fname_str(fname);
//         message_abort(cerr, "ERROR: Cannot write file "+fname_str);
//      }
//      return fp;
//   }
//}

int main(int argc,char **argv) {

   /*--- read arguments ---*/

   //int fflag=0,oflag=0,sflag=0,cflag=0;
   char *fvalue=NULL,*ovalue=NULL,*svalue=NULL;//,*cvalue=NULL;
   bool guess = false, assemble = false;
   char argu_key;

   opterr = 0;
   
   //while((argu_key = getopt(argc, argv,"f:o:s:c:")) != -1)
   while((argu_key = getopt(argc, argv,"f:o:s:ga")) != -1)
      switch(argu_key) {
         case 'f':
            fvalue = optarg;
            break;
         case 'o':
            ovalue = optarg;
            break;
         case 's':
            svalue = optarg;
            break;
         case 'g':
            guess = true;
            break;
         case 'a':
            assemble = true;
            break;
         case '?':
            if(optopt == 'f') 
               message_abort(cerr, "ERROR: Option -f requires an argument");
            else if(optopt == 'o') 
               message_abort(cerr, "ERROR: Option -o requires an argument");
            else if(optopt == 's') 
               message_abort(cerr, "ERROR: Option -s requires an argument");
            else message_abort(cerr, "ERROR: Unknown option -"+ string(1,optopt));
            break;
         default:
            message_abort(cerr, "ERROR: Unknown error");
      }
   //FILE *fp_f = open_rfile(fvalue);
   //FILE *fp_o = open_wfile(ovalue);
   //catch(logic_error& e) {
   //   cerr << e.what() << '\n';
   //   cout << "in catch\n";
   //   abort();
   //}
   //FILE *fp_c;
   char chainid;
   if(!fvalue) message_abort(cerr, "ERROR: Option -f not found");
   if(!ovalue) message_abort(cerr, "ERROR: Option -o not found");
   if(!svalue) {
      cerr << "WARNING: No segid provided, using W as default\n";
      chainid = 'W';
   }
   else {
      chainid = svalue[0];
      if(svalue[1]!='\0') { 
      cerr << "WARNING: Two or more chars were provided as segid" <<
         ", using the first char " << svalue[0] << endl;
      }
   }
   //if(cvalue) fp_c = open_wfile(cvalue);
   //else       fp_c = stdout;



   /*--- end of read arguments ---*/

   /*--- read pdb ---*/
   PDB pdb(fvalue);
   PDBDef defo, defh, defhyd;
   defo.pushBack(PDBField::chainid,chainid);
   defh = defo;
   defo.pushBack(PDBField::atomtype,"O");
   defh.pushBack(PDBField::atomtype,"H");
   defhyd = defo;
   defhyd.pushBack(PDBField::resname,"H3O");
   if(guess) cout << "will geuss pdb information\n";
   if(assemble) cout << "will assemble water atoms\n";
   pdb.reorderWater(guess, true, assemble, defo, defh, defhyd);
   pdb.write2file(ovalue);

   defhyd.print();

   //fprintf(fp_c, "%zu\n", center);
   //if(cvalue) fclose(fp_c);
}
