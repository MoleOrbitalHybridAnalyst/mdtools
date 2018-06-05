#include <iostream>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <unistd.h>
#include "pdb.h"
#include "utili.h"

#include "mpi.h"

using namespace std;
using namespace PDB_NS;

void print_help(ostream& s) {
   s << "Usage:\n";
   s << "\tgendensity -g grid.xyz -p pdb_dir"<< 
      " -d density_dir -r resolution -n nframes [-t transform]\n";
}

void message_abort(ostream& s, string str, bool help) {
   s << str << "\n";
   if(help) print_help(s);
   abort();
}

void message_abort(ostream& s, string str) {
   message_abort(s, str, true);
}

template <class T, class S>
double fsw3(const T& dist, const S& resols) {
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

// pos =  pos \cdot trans
template <class T, class S>
inline
void transform(T& pos, const S& trans) {
   auto pos_copy = pos;
   pos[0] = 
      pos_copy[0] * trans[0] + pos_copy[1] * trans[3] + pos_copy[2] * trans[6];
   pos[1] = 
      pos_copy[0] * trans[1] + pos_copy[1] * trans[4] + pos_copy[2] * trans[7];
   pos[2] = 
      pos_copy[0] * trans[2] + pos_copy[1] * trans[5] + pos_copy[2] * trans[8];
}

int main(int argc, char **argv) {
   char *gvalue = NULL, *pvalue = NULL, *dvalue = NULL;
   char *nvalue = NULL, *rvalue = NULL;
   char *tvalue = NULL;
   int  nframes;
   float reso;
   
   opterr = 0;
   char argu_key;
   while((argu_key = getopt(argc, argv,"g:p:d:n:r:t:")) != -1)
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
         case 't':
            tvalue = optarg;
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
            else if(optopt == 't') 
               message_abort(cerr, "ERROR: Option -t requires an argument");
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
   //vector<array<float, 3>> grids;
   vector<Vector> grids;
   ifstream fsgrid(gvalue);
   if(!fsgrid.is_open())
      message_abort(cerr, "ERROR: cannot open grid file", false);
   string line;
   while(getline(fsgrid, line)) {
      stringstream ss(line);
      //array<float, 3> grid;
      Vector grid(0.0);
      ss >> grid[0] >> grid[1] >> grid[2];
      grids.push_back(grid);
   }

   /*--- if tvalue is set, then read transformation info ---*/
   vector<vector<int> > transform_base;
   // transform matrices stored in row major order
   // atom position is considered to be a row vector
   // thus new_position = old_position \cdot matrix
   vector<array<double, 9> > transformations;
   // number of transformations will act on atoms
   unsigned ntransformations = 0;
   // number of transformation types
   unsigned ntranstypes = 0;
   if(tvalue) {
      /*--- read transform base ---*/
      ifstream fs_transbase((string(tvalue) + ".base").c_str());
      if(!fs_transbase.is_open())
         message_abort(cerr, "ERROR: cannot open transform base", false);
      string line;
      while(getline(fs_transbase, line)) {
         stringstream ss(line);
         int tbase;
         transform_base.push_back(vector<int>());
         while(ss >> tbase) {
            transform_base.back().push_back(tbase);
         }
         if(
               transform_base[0].size() !=
               transform_base.back().size()
           ) {
            cerr << "ERROR: inconsistent column numbers in "
               << tvalue << ".base";
            abort();
         }
      }
      ntranstypes = transform_base[0].size();
      ntransformations = transform_base.size();

      /*--- read transformations ---*/
      transformations.resize(ntranstypes);
      for(unsigned i = 0; i < ntranstypes; ++i) {
         ifstream fs_trans((string(tvalue) + to_string(i)).c_str());
         if(!fs_trans.is_open()) {
            string errmsg = string(tvalue) + to_string(i);
            errmsg = "ERROR: cannot open " + errmsg;
            message_abort(cerr, errmsg, false);
         }
         string line;
         unsigned index_in_matrix;
         double matrix_entry;
         while(getline(fs_trans, line)) {
            stringstream ss(line);
            while(ss >> matrix_entry) {
               transformations[i][index_in_matrix++] = matrix_entry;
               if(index_in_matrix > 9) 
                  cerr << "ERROR: 10th entry of a matrix is provided in "
                     << tvalue << i;
            }
         }
      }

   } // end of if(tvalue)

   // intialize mpi
   MPI_Init(NULL, NULL);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


   /*--- read pdbs and construct density ---*/
   for(int num = world_rank; num  < nframes; num += world_size) {

      stringstream sspdb;
      sspdb << string(pvalue) <<'/'<< num << ".pdb";
      cout << sspdb.str() << endl;

      PDB pdb(sspdb.str());
      
      /*--- construct density ---*/
      stringstream ssden;
      // will be the file name of the generated den
      unsigned index_den_start = num;
      unsigned index_den_end = num + 1;
      if(tvalue) {
         index_den_start *= ntransformations;
         index_den_end *= ntransformations;
      }
      for(unsigned index_den = index_den_start; 
            index_den < index_den_end; ++ index_den) 
      {
         ssden << string(dvalue) <<'/'<< index_den <<".dat";
         ofstream fsden(ssden.str());
         if(!fsden.is_open()) {
            string errmsg = 
               string(dvalue) + '/' + to_string(index_den) + ".dat";
            errmsg = "ERROR: cannot open " + errmsg;
            message_abort(cerr, errmsg, false);
         }

         fsden << "#! FIELDS density" << endl;
         fsden << fixed;
         
         for(auto grid_iter = grids.begin(); 
               grid_iter != grids.end(); ++grid_iter)
         {
            double den = 0;
            for(size_t index = 0; index < pdb.getNatoms(); ++index) {

               // the position of atom index
               auto pos = pdb.getCoordinates(index);

               /*--- do the transformations if any ---*/
               if(tvalue) {
                  for(int num_trans : 
                        transform_base[index_den - index_den_start])
                     for(int i_trans = 0; i_trans < num_trans; ++i_trans) {
                        transform(pos, 
                              transformations[index_den - index_den_start]);
                     }
               }

               den += fsw3(
                     pdb.pbcDistance(*grid_iter, pos), resols);

            }
            fsden << ' ' << setprecision(8) << den << '\n';
         }

         fsden.close();

         
      } // end of loop in transformations
      
   } // end of mpi loop in pdbs

   MPI_Finalize();

   return 0;
   
}
