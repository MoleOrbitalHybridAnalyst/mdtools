# assume q mol x y z are in trj file
# unwrap coordinates using molid
from lammpstrj import dump, np
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('lammpstrj', help = 'input trj')
    parser.add_argument('dipole_out', help = 'output file')
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

    fp_out = open(args.dipole_out, 'w')
    print('#! FIELDS timestep dx dy dz', file = fp_out)

    # store file name but not read
    d = dump(args.lammpstrj, 0)

    # I don't know why Pizza.py used 'rb' when reading in next()
    # here I just manually do partial work in read_all()
    # and compute dipole on the fly of reading
    for fn in d.flist:
        if fn[-3:] == '.gz':
            f = popen("%s -c %s" % (gunzip,fn),'r')
        else: f = open(fn)

        snap =  d.read_snapshot(f)
        while snap:

            dipole = np.zeros(3)
            box = make_box(snap)

            #--- do dipole calculatin ---#
            molid2atomids = {}
            natoms = snap.natoms

            q = d.names['q']
            x = d.names['x']
            y = d.names['y']
            z = d.names['z']
            mol = d.names['mol']

            # construct molids2atomids
            for iatom in range(natoms):
                molid_this_atom = int(snap.atoms[iatom, mol])
                if molid_this_atom in molid2atomids:
                    molid2atomids[molid_this_atom].append(iatom)
                else:
                    molid2atomids[molid_this_atom] = [iatom]

            # unwrap coords according to molids
            for molid, atomids in molid2atomids.items():
                id_ref = atomids[0]
                r_ref = snap.atoms[id_ref, [x, y, z]]
                dipole += snap.atoms[id_ref, q] * r_ref

                # @@@
                #if len(atomids) != 3:
                #    raise RuntimeError("broken water", snap.time, atomids)
                # @@@
                for atomid in atomids[1:]:
                    diff = snap.atoms[atomid, [x, y, z]] - r_ref
                    diff = pbc(diff, box)
                    # @@@
                    #if(sum(diff**2) > 4):
                    #    raise RuntimeError("large bond", snap.time, diff, box)
                    # @@@
                    dipole += snap.atoms[atomid, q] * (diff + r_ref)

            # print resulting dipole
            print(" %d %8f %8f %8f"\
                    %(snap.time, dipole[0], dipole[1], dipole[2]),
                    file = fp_out)

            snap = d.read_snapshot(f)

        f.close()

    fp_out.close()
