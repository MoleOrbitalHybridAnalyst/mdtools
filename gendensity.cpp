#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include <unistd.h>
#include "pdb.h"

#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace PDB_NS;

void print_help(ostream& s) {
   s << "Usage:\n";
   s << "\tgendensity -g grid.xyz -p pdb_dir"<< 
      " -d density_dir -r resolution -n nframes\n";
}

void message_abort(ostream& s, string str) {
   s << str << "\n";
   print_help(s);
   abort();
}

double fsw3(const array<float, 3>& dist, const array<float ,3>& resols) {
   double den = 1.0;
   for(int i = 0; i < 3; ++i) {
      if(dist[i] > resols[i] or dist[i] < -resols[i]) {
         return 0.0;
      } else {
         double t = fabs(dist[i] / resols[i]) -0.5;
         //den *= (t * (t*t*(-6*t*t+5)-1.875) + 0.5);
         den *= (t*(2*t*t - 1.5) + 0.5);
      }
   }
   return den;
}

int main(int argc, char **argv) {
   char *gvalue=NULL, *pvalue=NULL, *dvalue=NULL;
   char *nvalue=NULL, *rvalue=NULL;
   int  nframes;
   float reso;
   
   opterr = 0;
   char argu_key;
   while((argu_key = getopt(argc, argv,"g:p:d:n:r:")) != -1)
      switch(argu_key) {
         case 'g':
            gvalue = optarg;
            break;
         case 'p':
            pvalue = optarg;
            break;
         case 'd':
            dvalue = optarg;
            break;
         case 'n':
            nvalue = optarg;
            break;
         case 'r':
            rvalue = optarg;
            break;
         case '?':
            if(optopt == 'g') 
               message_abort(cerr, "ERROR: Option -g requires an argument");
            else if(optopt == 'p') 
               message_abort(cerr, "ERROR: Option -p requires an argument");
            else if(optopt == 'd') 
               message_abort(cerr, "ERROR: Option -d requires an argument");
            else if(optopt == 'n') 
               message_abort(cerr, "ERROR: Option -n requires an argument");
            else if(optopt == 'r') 
               message_abort(cerr, "ERROR: Option -r requires an argument");
            else message_abort(cerr, "ERROR: Unknown option -"+ string(1,optopt));
            break;
         default:
            message_abort(cerr, "ERROR: Unknown error");
      }

   // read grid
   //float resox, resoy ,resoz;
   //resox = resoy = resoz = reso;
   if(!gvalue) {
      message_abort(cerr, "ERROR: Option -g not found");
   }
   if(!pvalue) {
      message_abort(cerr, "ERROR: Option -p not found");
   }
   if(!dvalue) {
      message_abort(cerr, "ERROR: Option -d not found");
   }
   if(!nvalue) {
      message_abort(cerr, "ERROR: Option -n not found");
   }
   if(!rvalue) {
      message_abort(cerr, "ERROR: Option -r not found");
   }
   stringstream(nvalue) >> nframes;
   stringstream(rvalue) >> reso;
   array<float, 3> resols;
   resols[0] = reso; resols[1] = reso; resols[2] = reso;
   vector<array<float, 3>> grids;
   ifstream fsgrid(gvalue);
   string line;
   while(getline(fsgrid, line)) {
      stringstream ss(line);
      array<float, 3> grid;
      ss >> grid[0] >> grid[1] >> grid[2];
      grids.push_back(grid);
   }

   // read pdbs and construct density
   for(int num = 0; num  < nframes; ++num) {
      stringstream sspdb;
      sspdb << string(pvalue) <<'/'<< num << ".pdb";
      cout << sspdb.str() << endl;

      //auto t1 = high_resolution_clock::now();
      PDB pdb(sspdb.str());
      //auto t2 = high_resolution_clock::now();
      //cout << duration_cast<duration<double>>(t2-t1).count() << endl;
      // construct density
      stringstream ssden;
      ssden << string(dvalue) <<'/'<< num <<".dat";
      ofstream fsden(ssden.str());
      //t1 = high_resolution_clock::now();
      for(auto grid_iter = grids.begin(); grid_iter != grids.end(); ++grid_iter)
      {
         double den = 0;
         for(size_t index = 0; index < pdb.getNatoms(); ++index) {
            //array<float, 3> coord;
            //array<float, 3> coord = pdb.getCoordinates(index);
            //coord[0] = pdb.getX(index);
            //coord[1] = pdb.getY(index);
            //coord[2] = pdb.getZ(index);
//            coord[0] = pdb.getX(index) - 79.63125133514404;
//            coord[1] = pdb.getY(index) - 42.08624982833863;
//            coord[2] = pdb.getZ(index) - 112.31375217437744;
            //den += fsw3(pdb.pbcDistance(*grid_iter, coord), resols);
            den += fsw3(
                  pdb.pbcDistance(*grid_iter, pdb.getCoordinates(index)),
                  resols);

//            if(grid_iter==grids.begin()+467) {
//               //cout << coord[0] << endl;
//               if(index==1) {
//                  printf("%f %f %f\n%f %f %f\n",
//                        coord[0], coord[1], coord[2],
//                        (*grid_iter)[0],(*grid_iter)[1],(*grid_iter)[2]);
//               }
//               cout << den << endl;
//            }
         }
         fsden << den << '\n';
      }
      fsden.close();
      //t2 = high_resolution_clock::now();
      //cout << duration_cast<duration<double>>(t2-t1).count() << endl;
   }
   
}
