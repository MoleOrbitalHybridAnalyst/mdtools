#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

//#include <experimental/array>

#include <unistd.h>
#include "pdb.h"
#include "utili.h"

#include "mpi.h"

#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace PDB_NS;

void print_help(ostream& s) {
   s << "Usage:\n";
   s << "\tgendensity -g grid.xyz -p pdb_dir"<< 
      " -d density_dir -r resolution -n nframes -c\n";
}

void message_abort(ostream& s, string str) {
   s << str << "\n";
   print_help(s);
   abort();
}

template <class T, class S>
double fsw3(const T& dist, const S& resols) {
   double den = 1.0;
   for(int i = 0; i < 3; ++i) {
      if(dist[i] > resols[i] or dist[i] < -resols[i]) {
         return 0.0;
      } else {
         double t = fabs(dist[i] / resols[i]) -0.5;
         den *= (t * (t*t*(-6*t*t+5)-1.875) + 0.5);
      }
   }
   
   return den;
}

template <class T, class S, class Q>
double fsw3(const T& dist, const S& resols, Q& derivs) {
   Vector f, df;
   for(int i = 0; i < 3; ++i) {
      if(dist[i] > resols[i] or dist[i] < -resols[i]) {
         return 0.0;
      } else {
         double t = fabs(dist[i] / resols[i]) -0.5;
         //den *= (t * (t*t*(-6*t*t+5)-1.875) + 0.5);
         f[i] = (t*(2*t*t - 1.5) + 0.5);
         df[i] = 
            static_cast<int>((dist[i] > 0.0) - (dist[i] < 0.0)) * (6*t*t - 1.5);
      }
   }
   
   derivs[0] = df[0] *  f[1] *  f[2];
   derivs[1] =  f[0] * df[1] *  f[2];
   derivs[2] =  f[0] *  f[1] * df[2];

   return f[0] * f[1] * f[2];
}

int main(int argc, char **argv) {
   char *gvalue=NULL, *pvalue=NULL, *dvalue=NULL;
   char *nvalue=NULL, *rvalue=NULL;
   int  nframes;
   float reso;
   bool do_derivs = false;
   
   opterr = 0;
   char argu_key;
   while((argu_key = getopt(argc, argv,"g:p:d:n:r:c")) != -1)
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
         case 'c':
            do_derivs = true;
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
   if(do_derivs) 
      cout << "will calculate derivs\n";
   else 
      cout << "will not calculate derivs\n";
   stringstream(nvalue) >> nframes;
   stringstream(rvalue) >> reso;
   array<float, 3> resols;
   resols[0] = reso; resols[1] = reso; resols[2] = reso;
   //vector<array<float, 3>> grids;
   vector<Vector> grids;
   ifstream fsgrid(gvalue);
   string line;
   while(getline(fsgrid, line)) {
      stringstream ss(line);
      //array<float, 3> grid;
      Vector grid(0.0);
      ss >> grid[0] >> grid[1] >> grid[2];
      grids.push_back(grid);
   }

   // intialize mpi
   MPI_Init(NULL, NULL);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


   // read pdbs and construct density
   for(int num = world_rank; num  < nframes; num += world_size) {
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
      //write the header of COLVAR
      fsden << "#! FIELDS density" << endl;
      fsden << fixed;
      //t1 = high_resolution_clock::now();
      for(auto grid_iter = grids.begin(); grid_iter != grids.end(); ++grid_iter)
      {
         double den = 0;
         for(size_t index = 0; index < pdb.getNatoms(); ++index) {

            Vector derivs;
            if(do_derivs) 
               den += fsw3(
                     pdb.pbcDistance(*grid_iter, pdb.getCoordinates(index)),
                     resols, derivs);
            else
               den += fsw3(
                     pdb.pbcDistance(*grid_iter, pdb.getCoordinates(index)),
                     resols);

         }
         fsden << ' ' << setprecision(8) << den << '\n';
      }
      fsden.close();
      
   }

   // finialize MPI
   MPI_Finalize();
   
}
