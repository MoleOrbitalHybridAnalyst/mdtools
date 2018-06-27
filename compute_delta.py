# read evb.out and lammpstrj
# generate { [timestep, pivot_oxygen_id, i, c_i/c_1, delta_i, r_O_1,O_i] }

# Assumptions:
# (1) id, mol and type are provided in trj
# (2) ids are consecutive numbers starting from 1
# (3) first solvation shell oxygen always accepts the closest hydronium hydrogen

# Behaviour expected if assumptions are violated:
# (1) Crash due to ValueError
# (2) Crash due to IndexError
# (3) Unmatched i, delta_i

# Considerations:
# (1) 'pivot' means the state with largest eig; not the search center

# Behaviour expected if considerations are violated:
# (1) 'pivot' can have only two hydrogens resulting in very wrong c_i/c_1,
#     delta_i pair because the missing hydrogen must be the bridge hydrogen
#     for some acceptor oxygen

from evb_out import *
from lammpstrj import *
import argparse 

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('evb_out_fn', help = 'input evb.out')
    parser.add_argument('lammpstrj', help = 'input lammps trj')
    parser.add_argument('output_dat', help = 'output results')
    parser.add_argument('--o_types', 
            help = 'lammps types for oxygens',
            default = '3,5' )
    return parser.parse_args()

def pbc(x, box):
    n = (x + 0.5 * box) // box
    return x - n * box

def pbc_distance2(iatom, jatom, trj, box):
    ri = np.array( trj.atom(iatom, 'x', 'y', 'z') ).flatten()
    rj = np.array( trj.atom(jatom, 'x', 'y', 'z') ).flatten()
    return np.sum( pbc(rj - ri, box) **2 )


if __name__ == '__main__':

    #--- parse argus ---#
    args = parse()
    fpout = open(args.output_dat, 'w')
    print('#! FIELDS timestep donor_atom_serial acceptor_atom_serial' + 
            ' ci_ratio delta r_OO', file = fpout)
    o_types = np.vectorize(int)(args.o_types.split(',')).tolist()
    
    #--- read evb.out and lammpstrj ---#
    evb = evb_out(args.evb_out_fn)
    print("finish reading", args.evb_out_fn)
    trj = dump(args.lammpstrj)

    trj_timesteps = trj.time()
    evb_timesteps = evb.time_trj.tolist()

    for timestep in trj_timesteps:

        trj.tselect.one(timestep)

        #--- get box, natoms in this frame ---#
        natoms = trj.snaps[0].natoms
        box = np.zeros(3)
        box[0] = trj.snaps[0].xhi - trj.snaps[0].xlo
        box[1] = trj.snaps[0].yhi - trj.snaps[0].ylo
        box[2] = trj.snaps[0].zhi - trj.snaps[0].zlo


        try:
            indx_in_evb_out = evb_timesteps.index(timestep)
        except ValueError:
            raise Exception(
                    "timestep %d is not found in"%timestep,
                    args.evb_out_fn )

        evb_state = evb.states_trj[indx_in_evb_out]

        #-- get the evb coeffs --#
        eig = evb.eigenv_trj[indx_in_evb_out]
        coeff_pivot = max(eig)
        pivot_molid = evb_state.mol_B[np.argmax(eig)]
        if np.argmax(eig):
            print('timestep %d worth being checked'%timestep)
        coeffs_first_solvshell = []

        #-- find the molids to be used --#
        first_solvshell_molids = []
        for coeff, molid_A, molid_B in zip(eig, \
                evb_state.mol_A, evb_state.mol_B):
            # if molid_B's parent is pivot_molid, then it's in first solvshell
            if molid_A == pivot_molid:
                first_solvshell_molids.append(molid_B)
                coeffs_first_solvshell.append(coeff)
            # if molid_A's child is pivot_molid, then it's in first solvshell
            # if molid_A == -1 then molid_B is pivot not 1st solvshell
            if molid_B == pivot_molid and molid_A > 0:
                first_solvshell_molids.append(molid_A)
                # find the coeff of state of molid_A
                indx = evb_state.mol_B.tolist().index(molid_A)
                coeffs_first_solvshell.append(eig[indx])

        #-- get the atomids to be used --#
        molid2atomid = {}
        for iatom in range(1, natoms + 1):

            molid_of_iatom = int(trj.atom(iatom, 'mol')[0])

            # test if molid_of_iatom is in first solvshell
            try:
                first_solvshell_molids.index(molid_of_iatom)
            except ValueError:
                # is molid_of_iatom the pivot molid ?
                if molid_of_iatom != pivot_molid: continue

            if molid_of_iatom in molid2atomid:
                molid2atomid[molid_of_iatom].append(iatom)
            else:
                molid2atomid[molid_of_iatom] = [iatom]

        #-- find the oxygen and hydrogens in hydronium --#
        pivot_oxygen_id = -1
        pivot_hydrogen_ids = []
        for iatom in molid2atomid[pivot_molid]:
            # check if this is oxygen
            try:
                type_of_iatom = int(trj.atom(iatom, 'type')[0]) 
                o_types.index(type_of_iatom)
                pivot_oxygen_id = iatom
            except ValueError:
                pivot_hydrogen_ids.append(iatom)
        if pivot_oxygen_id == -1:
            raise Exception('does not find oxygen in hydronium' +
                    '(molid = %d) in step %d'%(pivot_molid, timestep))
        if len(pivot_hydrogen_ids) != 3:
            raise Exception('number of hydrogens in hydronium' +
                    '(molid = %d) is %d != 3 in step %d'\
                    %(pivot_molid, len(pivot_hydrogen_ids), timestep))

        #-- do the real calculation --#
        for coeff, molid in zip(coeffs_first_solvshell, first_solvshell_molids):
            # find the oxygen in the mol
            acceptor_oxygen_id = -1
            for iatom in molid2atomid[molid]:
                try:
                    type_of_iatom = int(trj.atom(iatom, 'type')[0])
                    o_types.index(type_of_iatom)
                    acceptor_oxygen_id = iatom
                except ValueError:
                    continue
            if acceptor_oxygen_id == -1:
                raise Exception('does not find oxygen in water' +
                     ' %d in step %d'%(molid, timestep))
            r_OO = np.sqrt(pbc_distance2(
                pivot_oxygen_id, acceptor_oxygen_id, trj, box) )
            # @@@@
            #print(trj.atom(pivot_oxygen_id, 'x', 'y', 'z'))
            #print(trj.atom(acceptor_oxygen_id, 'x', 'y', 'z'))
            #print(box)
            #print(r_OO)

            # find the closest hydronium hydrogen to the acceptor oxygen
            r_O2H = np.sum(box**2)
            bridge_hydrogen_id = -1
            for iatom in pivot_hydrogen_ids:
                dist2 = pbc_distance2(iatom, acceptor_oxygen_id, trj, box)
                if dist2 <= r_O2H:
                    bridge_hydrogen_id = iatom
                    r_O2H = dist2
            if bridge_hydrogen_id == -1:
                raise Exception('does not find bridge hydrogen for water' +
                        ' %d in step %d'%(molid, timestep))
            r_O2H = np.sqrt(r_O2H)
            
            # compute delta
            r_O1H = pbc_distance2(pivot_oxygen_id, bridge_hydrogen_id, trj, box)
            r_O1H = np.sqrt(r_O1H)
            delta = r_O2H - r_O1H

            # print out results for this mol
            print("", timestep, pivot_oxygen_id, acceptor_oxygen_id, 
                    coeff / coeff_pivot, delta, r_OO, file = fpout)


    fpout.close()
