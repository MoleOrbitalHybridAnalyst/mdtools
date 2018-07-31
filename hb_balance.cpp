// hb_balance trj r_oo theta
// will generate delta.dat and hb.dat
#include "lammpstrj.h"

using namespace LAMMPSTRJ_NS;

int main(int argc, char *argv[1]) {

   (void) argc;

   LAMMPSTRJ lmptrj(argv[1]);

   // @@@@
   std::ofstream fs_delta("delta.dat");
   std::ofstream fs_hb("hb.dat");
   fs_delta << 
      "#! FIELDS delta1 delta2 delta3 OTH HTH1 HTH2 HTH3 OW1 OW2 OW3\n";
   fs_hb << "#! FIELDS hb_out1 hb_in1 hb_out2 hb_in2 hb_out3 hb_in3\n";

   double r_oo, theta;
   std::string dstr(argv[2]);
   std::stringstream(dstr) >> r_oo;
   dstr = std::string(argv[3]);
   std::stringstream(dstr) >> theta;

   constexpr unsigned type_water_o = 1, type_water_h = 2, 
                      type_hyd_o = 3, type_hyd_h = 4;

   constexpr int max_ip = 3; // hydronium has 3 protons

   for(unsigned iframe = 0; iframe < lmptrj.size(); ++iframe) {

      fs_delta << ' ';

      // find the OTH HTH all OWs
      unsigned oth_index;
      std::vector<unsigned> all_ow_indexes, hth_indexes, all_hw_indexes;

      auto& snap = lmptrj[iframe];
      for(unsigned iatom = 0; iatom < snap.get_natoms(); ++iatom) {
         unsigned this_type = snap.get_type(iatom);
         if(this_type == type_water_h) all_hw_indexes.push_back(iatom);
         else if(this_type == type_water_o) all_ow_indexes.push_back(iatom);
         else if(this_type == type_hyd_h) hth_indexes.push_back(iatom);
         else if(this_type == type_hyd_o) oth_index = iatom;
      }

      std::vector<unsigned> true_ow_indexes;
      for(unsigned iatom_hth : hth_indexes) {

         // find the associated iatom_ow for iatom_hth
         double mindist2 = std::numeric_limits<double>::max();
         unsigned true_ow_index = all_ow_indexes[0];
         for(unsigned iatom_ow : all_ow_indexes) {
            double dist2 = snap.distance2(iatom_hth, iatom_ow);
            if(dist2 < mindist2) {
               mindist2 = dist2; true_ow_index = iatom_ow;
            }
         }
         true_ow_indexes.push_back(true_ow_index);

         // compute deltas
         fs_delta << 
            fabs(std::sqrt(mindist2) - snap.distance(iatom_hth, oth_index))
                  << ' ';
      }
      fs_delta << snap.get_id(oth_index) << ' '
               << snap.get_id(hth_indexes[0]) << ' '
               << snap.get_id(hth_indexes[1]) << ' '
               << snap.get_id(hth_indexes[2]) << ' '
               << snap.get_id(true_ow_indexes[0]) << ' '
               << snap.get_id(true_ow_indexes[1]) << ' '
               << snap.get_id(true_ow_indexes[2]) << '\n';
         
      // get the molids of true OWs
      std::vector<unsigned> true_ow_mols;
      for(unsigned iatom : true_ow_indexes)
         true_ow_mols.push_back(snap.get_mol(iatom));

      // find the bonding HW for true OWs
      std::vector<unsigned> true_hw_indexes[max_ip];
      for(unsigned iatom : all_hw_indexes) {
         unsigned this_mol_hw = snap.get_mol(iatom);
         if(this_mol_hw == true_ow_mols[0])
            true_hw_indexes[0].push_back(iatom);
         else if(this_mol_hw == true_ow_mols[1])
            true_hw_indexes[1].push_back(iatom);
         else if(this_mol_hw == true_ow_mols[2])
            true_hw_indexes[2].push_back(iatom);
      }

      //printf("OWs: %u %u %u HWs: %d %d %d %d %d %d\n",
      //      true_ow_indexes[0], true_ow_indexes[1], true_ow_indexes[2],
      //      true_hw_indexes[0][0], true_hw_indexes[1][0], true_hw_indexes[2][0],
      //      true_hw_indexes[0][1], true_hw_indexes[1][1], true_hw_indexes[2][1]);

      unsigned num_hb_out[max_ip];
      unsigned num_hb_in[max_ip];
      double r_oo2 = r_oo * r_oo;
      double cos_theta = cos(theta / 180 * 3.1415926535);
      for(int ip = 0; ip < max_ip; ip++) {

         num_hb_out[ip] = 0;
         num_hb_in[ip] = 0;

         unsigned first_shell_ow_index = true_ow_indexes[ip];

         // check whether OTH HTH true_OW is a HB in
         if(
            snap.distance2(oth_index, first_shell_ow_index) <= r_oo2 &&
            snap.cos_angle(hth_indexes[ip], oth_index, first_shell_ow_index)
            >= cos_theta
           )  num_hb_in[ip] ++;

         std::vector<unsigned> second_ows;
         std::vector<unsigned> second_wat_mols;

         // get OWs within r_oo of true_OWs and get HBs out
         for(unsigned ow_index : all_ow_indexes) {

            if(snap.distance2(ow_index, first_shell_ow_index) <= r_oo2
                  && ow_index != first_shell_ow_index) {

               second_ows.push_back(ow_index);
               second_wat_mols.push_back(snap.get_mol(ow_index));

               for(unsigned first_shell_hw_index : true_hw_indexes[ip]) {

                  // compute cos(angle) 
                  if(snap.cos_angle(first_shell_hw_index,
                                    first_shell_ow_index, ow_index) >= cos_theta
                    )  {
//                     printf("%u-%u...%u\n", 
//                          first_shell_ow_index, first_shell_hw_index, ow_index);
                     num_hb_out[ip] ++;
                  }
               } // loop for first shell water HWs
            }

         } // loop for all OWs


         // get second HWs from second_wat_mols
         std::vector<std::vector<unsigned>> second_hws(second_ows.size());
         for(unsigned iatom = 0; iatom < snap.get_natoms(); ++iatom) {
            if(snap.get_type(iatom) == type_water_h) {

               unsigned this_mol = snap.get_mol(iatom);
               for(unsigned imol = 0; imol < second_wat_mols.size(); ++imol) {
                  if(this_mol == second_wat_mols[imol])
                     second_hws[imol].push_back(iatom);
               }

            }
         }

         // check whether second shell water donor HBs to first water
         for(unsigned imol = 0; imol < second_wat_mols.size(); ++imol) {
//            printf("OW = %u HW = %u %u\n",
//                  second_ows[imol], second_hws[imol][0], second_hws[imol][1]);
            // NOTE: second shell water already meets r_oo criteria
            for(unsigned second_hw_index : second_hws[imol]) {
               if(snap.cos_angle(
                  second_hw_index, second_ows[imol], 
                  first_shell_ow_index) >= cos_theta) {
                  num_hb_in[ip] ++;
//                  printf("%u-%u...%u\n", 
//                       second_ows[imol], second_hw_index, first_shell_ow_index);
               }
            }
         }

      } // loop for ip

      fs_hb << ' ' <<
            num_hb_out[0] << ' ' << num_hb_in[0] << ' ' <<
            num_hb_out[1] << ' ' << num_hb_in[1] << ' ' <<
            num_hb_out[2] << ' ' << num_hb_in[2] << '\n';

   }

   // @@@@
   fs_delta.close();
   fs_hb.close();

   return 0;
}
