# assume q mol x y z are in trj file
# unwrap coordinates using molid

# Two things I am uncomfortable with in this code:
# (1) Coordinates are unwrapped by molid; this may be problematic with EVB
# (2) Bridge hydrogen are determined by looking for the closest hydrogen
from lammpstrj import dump, np
from evb_out import evb_out
import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('evb_out', help = 'input evb.out')
    parser.add_argument('lammpstrj', help = 'input trj(s)')
    parser.add_argument('dipole_out', help = 'output file')
    parser.add_argument('--q_water_o', help = 'water oxygen charge',
            default = -0.82)
    parser.add_argument('--q_water_h', help = 'water hydrogen charge',
            default = 0.41)
    parser.add_argument('--q_hyd_o', help = 'hydronium oxygen charge',
            default = -0.32)
    parser.add_argument('--q_hyd_h', help = 'hydronium oxygen charge',
            default = 0.44)
    parser.add_argument('--type_o', help = 'water/hydronium oxygen type',
            default = '1,3')
    parser.add_argument('--type_h', help = 'water/hydronium hydrogen type',
            default = '2,4')
    return parser.parse_args()

def pbc(d, box):
    n = (d + 0.5 * box) // box
    return d - n * box

def make_box(snap):
    box = []
    box.append(snap.xhi - snap.xlo)
    box.append(snap.yhi - snap.ylo)
    box.append(snap.zhi - snap.zlo)
    return np.array(box)

if __name__=='__main__':
    args = parse()

    q_water_o = float(args.q_water_o)
    q_water_h = float(args.q_water_h)
    q_hyd_o = float(args.q_hyd_o)
    q_hyd_h = float(args.q_hyd_h)

    type_o = np.vectorize(int)(args.type_o.split(',')).tolist()
    type_h = np.vectorize(int)(args.type_h.split(',')).tolist()

    fp_out = open(args.dipole_out, 'w')
    print('#! FIELDS timestep dx dy dz', file = fp_out)

    # store file name but not read
    d = dump(args.lammpstrj, 0)

    # read evb.out
    evb_out_fn = args.evb_out
    evb = evb_out(evb_out_fn)
    evb_time_trj = evb.time_trj.tolist()
    print("evb_out read")

    # I don't know why Pizza.py used 'rb' when reading in next()
    # here I just manually do partial work in read_all()
    # and compute dipole on the fly of reading
    for fn in d.flist:
        if fn[-3:] == '.gz':
            f = popen("%s -c %s" % (gunzip,fn),'r')
        else: f = open(fn)

        snap =  d.read_snapshot(f)
        while snap:

            print(snap.time)

            dipole = np.zeros(3)

            natoms = snap.natoms

            charges = np.zeros(natoms)
            weights = np.zeros(natoms)
            dipole = np.zeros(3)
            box = make_box(snap)

            q = d.names['q']
            x = d.names['x']
            y = d.names['y']
            z = d.names['z']
            mol = d.names['mol']
            typ = d.names['type'] 

            # construct mapping from molid to atomids
            molid2atomids = {}
            for iatom in range(natoms):
                molid_this_atom = int(snap.atoms[iatom, mol])
                if molid_this_atom in molid2atomids:
                    molid2atomids[molid_this_atom].append(iatom)
                else:
                    molid2atomids[molid_this_atom] = [iatom]

            # compute charges
            try:
                evb_timestep_indx = evb_time_trj.index(snap.time)
            except ValueError:
                raise Exception(
                        "cannot find timestep %d in evb_out"%snap.time)
            pivot_state_indx = np.argmax(evb.eigenv_trj[evb_timestep_indx])
            if pivot_state_indx != 0:
                # modify mol_A mol_B and eigenv so that evb.out 
                # is consistent with the trajectory
                evb.eigenv_trj[evb_timestep_indx][pivot_state_indx],\
                evb.eigenv_trj[evb_timestep_indx][0] = \
                evb.eigenv_trj[evb_timestep_indx][0],\
                evb.eigenv_trj[evb_timestep_indx][pivot_state_indx]

                evb.states_trj[evb_timestep_indx].mol_A[pivot_state_indx],\
                evb.states_trj[evb_timestep_indx].mol_A[0] = \
                evb.states_trj[evb_timestep_indx].mol_A[0],\
                evb.states_trj[evb_timestep_indx].mol_A[pivot_state_indx]

                evb.states_trj[evb_timestep_indx].mol_B[pivot_state_indx],\
                evb.states_trj[evb_timestep_indx].mol_B[0] = \
                evb.states_trj[evb_timestep_indx].mol_B[0],\
                evb.states_trj[evb_timestep_indx].mol_B[pivot_state_indx]

                evb.states_trj[evb_timestep_indx].mol_A[0] = -1
                evb.states_trj[evb_timestep_indx].mol_A[pivot_state_indx] = \
                evb.states_trj[evb_timestep_indx].mol_B[0]
                

            for state_indx, [mol_A, mol_B, eig] in enumerate(
                    zip(evb.states_trj[evb_timestep_indx].mol_A, 
                        evb.states_trj[evb_timestep_indx].mol_B, 
                        evb.eigenv_trj[evb_timestep_indx])
                    ):
                if mol_A == -1:
                    for iatom in molid2atomids[mol_B]:
                        try:
                            type_h.index(int(snap.atoms[iatom, typ]))
                            charges[iatom] += eig**2 * q_hyd_h
                            weights[iatom] += eig**2
                            # modify the charges of hydronium in trj
                            snap.atoms[iatom, q] = q_water_h
                        except ValueError:
                            try:
                                type_o.index(int(snap.atoms[iatom, typ]))
                                charges[iatom] += eig**2 * q_hyd_o
                                weights[iatom] += eig**2
                                # modify the charges of hydronium in trj
                                snap.atoms[iatom, q] = q_water_o
                            except ValueError:
                                raise Exception(
                                "wrong type of atom index %d in timestep %d"\
                                %(iatom, snap.time) )
                else:
                    # we need to find the third hydrogen for mol_B
                    # to do that, we first find the oxygen in mol_B
                    # then find the closest hydrogen in mol_A to that oxygen
                    oxygen_iatom = -1
                    for iatom in molid2atomids[mol_B]:
                        try:
                            type_h.index(int(snap.atoms[iatom, typ]))
                            charges[iatom] += eig**2 * q_hyd_h
                            weights[iatom] += eig**2
                        except ValueError:
                            try:
                                type_o.index(int(snap.atoms[iatom, typ]))
                                charges[iatom] += eig**2 * q_hyd_o
                                weights[iatom] += eig**2
                                oxygen_iatom = iatom
                            except ValueError:
                                raise Exception(
                                "wrong type of atom index %d in timestep %d"\
                                %(iatom, time) )
                    mindist2 = np.sum(box**2)
                    third_h_iatom = -1
                    for iatom in molid2atomids[mol_A]:
                        try:
                            type_h.index(int(snap.atoms[iatom, typ]))
                            dist2 = np.sum(
                                    pbc(
                                        snap.atoms[iatom, [x, y, z]] -
                                        snap.atoms[oxygen_iatom, [x, y, z]],
                                        box
                                            ) ** 2
                                        )
                            if dist2 < mindist2:
                                mindist2 = dist2
                                third_h_iatom = iatom
                        except ValueError:
                            pass

                    charges[third_h_iatom] += eig**2 * q_hyd_h
                    weights[third_h_iatom] += eig**2

            # finalize the charges
            total_charge = 0.0
            for iatom in range(natoms):
                charges[iatom] += snap.atoms[iatom, q] * (1.0 - weights[iatom])
                total_charge += charges[iatom]
            if abs(total_charge) > 1e-4:
                print("WARNING: total charge != 0 at timestep %d"%snap.time)

            # q * x = dipole
            for molid, atomids in molid2atomids.items():
                ref_iatom = atomids[0]
                for iatom in atomids:
                    try:
                        type_o.index(int(snap.atoms[iatom, typ]))
                        ref_iatom = iatom
                        break
                    except ValueError:
                        pass

                ref_x = snap.atoms[ref_iatom, [x, y, z]]
                for iatom in atomids:
                    pos = pbc(snap.atoms[iatom, [x, y, z]] - ref_x, box) + ref_x
                    dipole += charges[iatom] * pos

            print(" %d %.8f %.8f %.8f"\
                    %(snap.time, dipole[0], dipole[1], dipole[2]), 
                    file = fp_out)

            snap = d.read_snapshot(f)

    f.close()
