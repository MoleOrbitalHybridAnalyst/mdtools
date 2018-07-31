// special_pair trj delta.dat
// compute the life-time of non-transfer interval from the inital config
// count the number of new O_1x identities in the interval
// count the number of changes of O_1x identities in the interval
// sort the O_0-O_1x distance and output
#include "lammpstrj.h"

#include <fstream>
#include <set>
#include <algorithm>
#include <iostream>

using namespace LAMMPSTRJ_NS;


int main(int argc, char *argv[]) {

   (void) argc;

   LAMMPSTRJ lmptrj(argv[1], 0);
   std::fstream fs_delta(argv[2]);

   std::ofstream fs_spp("special_pair.dat");
   fs_spp << "#! FIELDS timestep num_new_o1x_ids num_changes_of_o1x"
          << " roo1x roo1y roo1z\n";

   std::string line;
   unsigned line_count = 0;
   unsigned pre_oth = static_cast<unsigned>(-1);
   std::set<unsigned> o1x_ids;
   unsigned num_changes = 0;
   unsigned pre_o1x = static_cast<unsigned>(-1);

   double d2ows[3];
   int sorted_indxs[3] = {0, 1, 2};

   while(std::getline(fs_delta, line)) {
      if(line[0] == '#') continue;
      line_count ++;
      std::stringstream ss(line);
      unsigned oth, ows[3];
      double junk;
      ss >> junk >> junk >> junk; ss >> oth;
      ss >> ows[0] >> ows[0] >> ows[0] >> ows[0] >> ows[1] >> ows[2];
      if(pre_oth == static_cast<unsigned>(-1)) pre_oth = oth;

      if(oth != pre_oth) break;

      /*--- find the O_1x using trj ---*/
      const auto& snap = lmptrj.read_snap();

      // find the index of oth and ows
      unsigned indx_oth, indx_ows[3];
      indx_oth = snap.find_atom_by_id(oth);
      indx_ows[0] = snap.find_atom_by_id(ows[0]);
      indx_ows[1] = snap.find_atom_by_id(ows[1]);
      indx_ows[2] = snap.find_atom_by_id(ows[2]);

      for(int ip = 0; ip < 3; ++ip)
         d2ows[ip] = snap.distance(indx_oth, indx_ows[ip]);
//      std::cout << d2ows[0] << ' ' << d2ows[1] << ' ' << d2ows[2] << '\n';
      //int min_pos = 
      //   static_cast<int>(std::min_element(d2ows, d2ows + 3) - d2ows);
      std::sort(sorted_indxs, sorted_indxs + 3, 
            [&d2ows](int i1, int i2){return d2ows[i1] < d2ows[i2];});
      unsigned o1x = ows[sorted_indxs[0]];
      o1x_ids.insert(o1x);
      if(pre_o1x == static_cast<unsigned>(-1)) pre_o1x = o1x;
      if(o1x != pre_o1x) num_changes ++;

      //std::cout << line_count << ' ' 
      fs_spp << ' ' << line_count << ' ' 
             << o1x_ids.size() << ' '
             << num_changes << ' '
             << d2ows[sorted_indxs[0]] << ' '
             << d2ows[sorted_indxs[1]] << ' '
             << d2ows[sorted_indxs[2]] << '\n';

      pre_oth = oth;
      pre_o1x = o1x;
   }

   //std::cout << line_count << ' ' 
   fs_spp << ' ' << line_count << ' ' 
          << o1x_ids.size() << ' '
          << num_changes << ' '
          << d2ows[sorted_indxs[0]] << ' '
          << d2ows[sorted_indxs[1]] << ' '
          << d2ows[sorted_indxs[2]] << '\n';

   fs_spp.close();
   fs_delta.close();

   return 0;
}
