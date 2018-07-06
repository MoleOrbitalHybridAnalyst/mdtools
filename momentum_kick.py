""" 
 assume mol type x y z vx vy vz are in input_lammpstrj
 assume type are consecutive integers starting from 1

 mode file format:
  id0 dx0 dy0 dz0
  id1 dx1 dy1 dz1
  id2 dx2 dy2 dz2

 note that using special atom ids are assuming there is only one proton
 behaviour of using special atom ids with multiple protons is not well-
 exceptionaled !

 special atom ids:
 OTH: atom whose type == type_hyd_o
 OW1: water oxygen which has strongest HB with OTH
 OW2: water oxygen which has second strongest HB with OTH
 OW3: water oxygen which has third strongest HB with OTH
 HTH1: bridge hydrogen between OTH and OW1
 HTH2: bridge hydrogen between OTH and OW2
 HTH3: bridge hydrogen between OTH and OW3
 HW1*: water hydrogen which is bonded to OW1
 HW2*: water hydrogen which is bonded to OW2
 HW3*: water hydrogen which is bonded to OW3

"""
# TODO 
# (1) smarter exception
# (2) multiple protons with atom selection


from lammpstrj import *
import argparse

import sys

class Parser:
    """ class of argument parser """

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('input_lammpstrj', help = 'input trj')
        parser.add_argument('mode', help = 'input mode file')
        parser.add_argument('output_lammpstrj', help = 'output trj')
        parser.add_argument('--timestep', '-t', 
                help = 'timestep to be selected in input trj',
                default = '-1')
        parser.add_argument('--coord_system', '-c', 
                help = 'timestep to be selected in input trj',
                default = 'x,y,z')
        parser.add_argument('--dT', 
                help = 'total temperature increase after kicking',
                default = '1.5')
        parser.add_argument('--hb_measure',
                help = 'descriptor of hydrogen bonds')
        parser.add_argument('--donot_unwrapping', '-u',
                help = 'do not unwrap when defining coordinate system',
                action = 'store_true' )
        parser.add_argument('--quiet', '-q',
                help = 'try to be quiet', action = 'store_true')
        parser.add_argument('--masses', '-m', required = True,
                help = 'list the atomic masses in the order of type')
        parser.add_argument('--mvv2e', help = 'energy converter',
                default = '2390.0573615334906')
        parser.add_argument('--boltz', help = 'boltzmann constant',
                default = '0.0019872067')
        parser.add_argument('--kick_type', 
                help = 'type of kick (preserve_direction or minimum_perturb)',
                default = 'preserve_direction')
        parser.add_argument('--type_water_o', help = 'type of water oxygen')
        parser.add_argument('--type_water_h', help = 'type of water hydrogen')
        parser.add_argument('--type_hyd_o', help = 'type of hydronium oxygen')
        parser.add_argument('--type_hyd_h', help = 'type of hydronium hydrogen')
        self.args = parser.parse_args()

class Parameters:
    """ hard coded parameters used in this code """

    def __init__(self, args):
        self.args = args

        # these should be got from update.cpp of lammps
        self.mvv2e = float(args.mvv2e)
        self.boltz = float(args.boltz)

        self.special_atomids = \
              [ 'OTH' , 'HTH1', 'HTH2', 'HTH3', \
                'OW1' , 'OW2' , 'OW3' ,         \
                'HW11', 'HW21', 'HW31',         \
                'HW12', 'HW22', 'HW32'  ]
        # atomids that require defined hb_measure
        self.hb_measure_atomids = \
              [ 'OW1' , 'OW2' , 'OW3' , \
                'HTH1', 'HTH2', 'HTH3', \
                'HW11', 'HW21', 'HW31', \
                'HW12', 'HW22', 'HW32'  ]
        # atomids that require defined type_water_o
        self.type_water_o_atomids = \
              [ 'OW1' , 'OW2' , 'OW3' , \
                'HTH1', 'HTH2', 'HTH3', \
                'HW11', 'HW21', 'HW31', \
                'HW12', 'HW22', 'HW32'  ]
        # atomids that require defined type_water_h
        self.type_water_h_atomids = \
              [ 'HW11', 'HW21', 'HW31', \
                'HW12', 'HW22', 'HW32'  ]
        # atomids that require defined type_hyd_o
        self.type_hyd_o_atomids = self.special_atomids
        # atomids that require defined type_hyd_h
        self.type_hyd_h_atomids = \
              [ 'OW1' , 'OW2' , 'OW3' , \
                'HTH1', 'HTH2', 'HTH3', \
                'HW11', 'HW21', 'HW31', \
                'HW12', 'HW22', 'HW32'  ]
        try: self.type_water_o = int(args.type_water_o)
        except TypeError: self.type_water_o = None
        try: self.type_water_h = int(args.type_water_h)
        except TypeError: self.type_water_h = None
        try: self.type_hyd_o = int(args.type_hyd_o)
        except TypeError: self.type_hyd_o = None
        try: self.type_hyd_h = int(args.type_hyd_h)
        except TypeError: self.type_hyd_h = None

    def translate_atomid(self, atomid, lmptrj):
        if type(atomid) is int: return atomid
        elif type(atomid) is str or type(atomid) is np.str_:

            global special_atomids_values
            if atomid in special_atomids_values:
                return special_atomids_values[atomid]

            if atomid not in self.special_atomids:
                raise RuntimeError('cannot understand atomid ' + atomid)

            # check required info was defined
            if atomid in self.hb_measure_atomids and \
                    hb_measure.measure_type is None:
                raise RuntimeError(
                    'using ' + str(atomid) + ' requires hb_measure')
            if atomid in self.type_water_o_atomids \
                    and self.type_water_o is None:
                raise RuntimeError(
                    'using ' + str(atomid) + ' requires type_water_o')
            if atomid in self.type_water_h_atomids \
                    and self.type_water_h is None:
                raise RuntimeError(
                    'using ' + str(atomid) + ' requires type_water_h')
            if atomid in self.type_hyd_o_atomids and self.type_hyd_o is None:
                raise RuntimeError(
                    'using ' + str(atomid) + ' requires type_hyd_o')
            if atomid in self.type_hyd_h_atomids and self.type_hyd_h is None:
                raise RuntimeError(
                    'using ' + str(atomid) + ' requires type_hyd_h')
    
            # real translation
            if atomid == 'OTH':
                # assuming only 1 OTH present
                iatom = lmptrj.get_iatom_from_type(self.type_hyd_o)
                if len(iatom) == 1:
                    special_atomids_values['OTH'] = iatom[0]
                else:
                    raise RuntimeError('mulitple OTH found')
                return iatom[0]
            elif atomid[:3] == 'HTH':
                # if associated OW(atomid[3]) is undefined, identify it
                # NOTE:
                # when identifying it, all HTH* and other two OW* 
                # are also identified
                atomid_ow = 'OW' + atomid[3]
                if atomid_ow not in special_atomids_values:
                    self.translate_atomid(atomid_ow, lmptrj)
                return special_atomids_values[atomid]
            elif atomid[:2] == 'OW':
                coord_oth = lmptrj.get_coord(
                        para.translate_atomid('OTH', lmptrj))

                #--- identify all OW* HTH* ---#
                # get bonded HTH to OTH
                iatoms_hth = lmptrj.get_iatom_from_type(self.type_hyd_h)
                hb_measures_hth = []
                iatoms_ow = []
                if len(iatoms_hth) != 3:
                    raise RuntimeError(\
                            'number of HTH is ' + \
                            str(len(iatoms_hth)) + ' != 3' )

                # find the OW*
                all_iatoms_ow = lmptrj.get_iatom_from_type(self.type_water_o)
                if len(all_iatoms_ow) < 3:
                    raise RuntimeError('number of all OW is less than 3')
                for iatom_hth in iatoms_hth:
                    coord_hth = lmptrj.get_coord(iatom_hth)
                            #para.translate_atomid('HTH', lmptrj))

                    ## find the best HB acceptor for iatom_hth
                    # find the closest HB acceptor for iatom_hth
                    iatom_best_ow = -1
                    mindist2 = np.inf
                    for iatom_ow in all_iatoms_ow:
                        #hbm = hb_measure.compute(
                        #        coord_oth, coord_hth, \
                        #        lmptrj.get_coord(iatom_ow))
                        #if hbm > best_hbm:
                        #    best_hbm = hbm; iatom_best_ow = iatom_ow
                        dist2 = pbc.distance2(coord_hth, lmptrj.get_coord(iatom_ow))
                        if dist2 < mindist2:
                            mindist2 = dist2; iatom_best_ow = iatom_ow
                    if iatom_best_ow == -1:
                        raise RuntimeError(
                                'cannot find best HB acceptor for bridge' +\
                                ' hydrogen ' + str(iatom_hth))
                    else:
                        iatoms_ow.append(iatom_best_ow)
                        best_hbm = hb_measure.compute(
                                coord_oth, coord_hth, \
                                lmptrj.get_coord(iatom_best_ow))
                        hb_measures_hth.append(best_hbm)

                # sort the iatoms_hth and iatoms_ow 
                # according to hb_measures_hth descendingly
                iatoms_ow = np.array(iatoms_ow)
                iatoms_hth = np.array(iatoms_hth)
                hb_measures_hth = np.array(hb_measures_hth)
                order = np.argsort(-hb_measures_hth)
                special_atomids_values['OW1'] = iatoms_ow[order][0]
                special_atomids_values['OW2'] = iatoms_ow[order][1]
                special_atomids_values['OW3'] = iatoms_ow[order][2]
                special_atomids_values['HTH1'] = iatoms_hth[order][0]
                special_atomids_values['HTH2'] = iatoms_hth[order][1]
                special_atomids_values['HTH3'] = iatoms_hth[order][2]

                return special_atomids_values[atomid]

            elif atomid[:2] == 'HW':
                iatom_ow = \
                    self.translate_atomid('OW' + atomid[2], lmptrj)
                iatoms_hw = \
                    lmptrj.get_iatom_from_molid(lmptrj.atom(iatom_ow, 'mol')[0])
                if len(iatoms_hw) != 2:
                    raise RuntimeError('number of HW is ' + \
                            str(len(iatoms_hw)) + ' != 2')
                special_atomids_values['HW' + atomid[2] + '1'] = iatoms_hw[0]
                special_atomids_values['HW' + atomid[2] + '2'] = iatoms_hw[1]

                return special_atomids_values[atomid]

        else:
            raise RuntimeError('int or str expected in translate_atomid')

class Box:
    """ class of simulation box """

    def __init__(self, snap):
        self.xhi = snap.xhi
        self.xlo = snap.xlo
        self.yhi = snap.yhi
        self.ylo = snap.ylo
        self.zhi = snap.zhi
        self.zlo = snap.zlo
        self.sizes = np.zeros(3)
        self.sizes[0] = self.xhi - self.xlo
        self.sizes[1] = self.yhi - self.ylo
        self.sizes[2] = self.zhi - self.zlo

class Pbc:
    """ class of periodic boundary condition """

    def __init__(self, snap):
        self.box = Box(snap)

    def shrink(self, x):
        """ shrink a vector to be within half box lengths """
        n = (x + 0.5 * self.box.sizes) // self.box.sizes
        return x - n * self.box.sizes

    def delta(self, x1, x2):
        """ compute a vector connecting from x1 to x2 with pbc respected """
        return self.shrink(x2 - x1)

    def distance2(self, x1, x2):
        """ compute the pbc distance square between x1 and x2 """
        return np.sum(self.delta(x1, x2)**2)

    def distance(self, x1, x2):
        """ compute the pbc distance between x1 and x2 """
        return np.sqrt(self.distance2(x1, x2))

    def unwrap(self, x, x_ref):
        """ unwrap x so that it is within half box lengths of x_ref """
        return x_ref + self.delta(x_ref, x)

class HB_Measure:

    def __init__(self, hb_measure_string):
        if hb_measure_string is None: 
            self.measure_type = None; return None
        measure_para = hb_measure_string.split(',')
        self.measure_type = measure_para[0]
        # parse measure_type
        if self.measure_type == 'delta':
            if len(measure_para) != 1:
                raise RuntimeError('no parameters needed for delta')
        elif measure_type == 'roo_delta':
            if len(measure_para) != 3:
                raise RuntimeError('roo and theta are needed for roo_delta')
            self.roo = float(measure_para[1])
            self.theta = float(measure_para[2])
        elif measure_type == 'roo_angle':
            raise RuntimeError('not implemented yet')
        elif measure_type == 'hb_balance':
            if len(measure_para) != 3:
                raise RuntimeError('roo and theta are needed for hb_balance')
            self.roo = float(measure_para[1])
            self.theta = float(measure_para[2])
        else:
            raise RuntimeError('cannot understand hb_measure', measure_type)

    def compute(self, x_o1, x_h, x_o2):
        """ compute the quanlity of HB; the larger the better """
        if self.measure_type == 'delta':
            return pbc.distance(x_o1, x_h) - pbc.distance(x_h, x_o2)
        elif self.measure_type == 'roo_delta':
            raise RuntimeError('not implemented yet')
        elif self.measure_type == 'roo_angle':
            raise RuntimeError('not implemented yet')
        elif self.measure_type == 'hb_balance':
            raise RuntimeError('not implemented yet')
        else:
            raise RuntimeError('cannot understand hb_measure', measure_type)



class CoordSystem:
    """ class of coordinate system """

    def __init__(self, coord_system_string, do_unwrapping = True):
        self.do_unwrapping = do_unwrapping
        self.ex = np.array([1.0, 0.0, 0.0])
        self.ey = np.array([0.0, 1.0, 0.0])
        self.ez = np.array([0.0, 0.0, 1.0])
        #--- read atom ids for coodinate system setup --#
        self.ids = coord_system_string.split(',')
        self.true_ids = []
        if len(self.ids) != 3:
            raise RuntimeError('3 atom ids expected for coordinate system')

    def construct(self):

        xyz_counter = 0
        for atomid in self.ids:
            # check if this is x, y, z
            if atomid == 'x' or atomid == 'y' or atomid == 'z':
                xyz_counter += 1
                continue
            self.true_ids.append(para.translate_atomid(atomid, d))

        if xyz_counter and xyz_counter != 3:
            raise RuntimeError(
                    'cannot mix x, y, z with atom ids in coord_system')
        if xyz_counter == 3: return

        self.x0 = d.get_coord(self.true_ids[0])
        self.x1 = d.get_coord(self.true_ids[1])
        self.x2 = d.get_coord(self.true_ids[2])

        # ex = id0->id1
        if self.do_unwrapping:
            self.ex = pbc.delta(self.x0, self.x1)
        else:
            self.ex = self.x1 - self.x0
        norm = np.sum(self.ex ** 2)
        if norm < 1e-100:
            raise RuntimeError('zero ex got')
        else:
            self.ex /= np.sqrt(norm)

        # ez = id0->id1 x id0->id2
        if self.do_unwrapping:
            self.ez = np.cross(self.ex, pbc.delta(self.x0, self.x2))
        else:
            self.ez = np.cross(self.ex, self.x2 - self.x0)
        norm = np.sum(self.ez ** 2)
        # if ez === \vec{0}, ez = (0, 0, 1.0) - ((0, 0, 1).ex) ex
        if norm < 1e-100:
            print('WARNING: zero ez got, use Gram-Schmidt instead', 
                    file = sys.stderr)
            self.ez = np.array([0,0,1]) - np.dot([0,0,1], self.ex) * self.ex
            norm = np.sum(self.ez ** 2)
        else:
            self.ez /= np.sqrt(norm)

        # ey = ez x ex
        self.ey = np.cross(self.ez, self.ex)

    def show(self):
        if len(self.true_ids):
            print(
             "Atom serials used in coordinate system definition: %d %d %d"\
                    %(self.true_ids[0], self.true_ids[1], self.true_ids[2]))
        print("Using coordinate System:")
        print("  ex = %.8f %.8f %.8f"\
                %(self.ex[0], self.ex[1], self.ex[2]))
        print("  ey = %.8f %.8f %.8f"\
                %(self.ey[0], self.ey[1], self.ey[2]))
        print("  ez = %.8f %.8f %.8f"\
                %(self.ez[0], self.ez[1], self.ez[2]))

class Mode:
    """ class of motion mode """

    def __init__(self, mode_string):
        data = np.genfromtxt(mode_string, dtype = str)
        self.ids = [line[0] for line in data]
        # movements along coordinate system
        self.movements = [np.vectorize(float)(line[1:4]) for line in data]

        # translate ids
        self.true_ids = []
        for atomid in self.ids:
            self.true_ids.append(para.translate_atomid(atomid, d))

        self.true_movements = None

    def compute(self, coord_system, dT, kick_type = 'preserve_direction'):
        # convert dT into dKinEng in mvv unit
        dKinEng = dof * dT * para.boltz / 2 / para.mvv2e

        # decompose movements onto true_movements
        self.true_movements = []
        for movement in self.movements:
            self.true_movements.append(
                 np.dot(
                     movement, 
                     [coord_system.ex, coord_system.ey, coord_system.ez]
                     )
                 )

        # normalize movements so that dT is achieved
        A = 0.0; B = 0.0
        for iatom, movement in zip(self.true_ids, self.true_movements):
            typ, vx, vy, vz = \
               np.array(d.atom(iatom, 'type', 'vx', 'vy', 'vz')).flatten()
            A += 0.5 * masses[int(typ)] * np.sum(movement**2)
            B += masses[int(typ)] * np.dot([vx, vy, vz], movement)
        #    A g^2 + B g == dKinEng
        # => g == (-B \pm √(B^2 + 4AdKinEng)) / (2A)
        g1 = (-B + np.sqrt(B**2 + 4 * A * dKinEng)) / (2 * A)
        g2 = (-B - np.sqrt(B**2 + 4 * A * dKinEng)) / (2 * A)
        # g1 is always positive so that movement direction conserved
        if kick_type == 'preserve_direction': g = g1
        elif kick_type == 'minimum_perturb':
            if abs(g1) < abs(g2): g = g1
            else: g = g2
        else: RuntimeError('cannot understand kick_type ' + kick_type)
        #print("dKinEng =",dKinEng,"g =", g, "B/A =", B/A)
        self.true_movements = np.array(self.true_movements) * g

    def show(self):
        if self.true_movements is None:
            print("Velocity change of each atom in coordinate system:")
            for iatom, movement in zip(self.true_ids, self.movements):
                print("  %d %8f %8f %8f"\
                  %(iatom, movement[0], movement[1], movement[2]))
        else:
            print("Velocity change of each atom in xyz coordinate system:")
            for iatom, movement in zip(self.true_ids, self.true_movements):
                print("  %d %f %f %f"\
                  %(iatom, movement[0], movement[1], movement[2]))


def get_coord(lmptrj, iatom):
    if len(lmptrj.time()) != 1:
        raise RuntimeError(
                "cannot use get_coord when multiple timesteps selected")
    return np.array(lmptrj.atom(iatom, 'x', 'y', 'z')).flatten()

def get_iatom_from_type(lmptrj, typ):
    """ gives atom ids in lmptrj whose type is typ """

    iatoms = []

    try: typ = int(typ)
    except TypeError:
        raise RuntimeError('cannot understand atom type ' + str(typ))
    for iatom in range(1, natoms + 1):
        if int(lmptrj.atom(iatom, 'type')[0]) == typ:
            iatoms.append(iatom)

    return iatoms

def get_iatom_from_molid(lmptrj, mol):
    """ gives atom ids in lmptrj whose type is typ """

    iatoms = []

    try: mol = int(mol)
    except TypeError:
        raise RuntimeError('cannot understand molid ' + str(typ))
    for iatom in range(1, natoms + 1):
        if int(lmptrj.atom(iatom, 'mol')[0]) == mol:
            iatoms.append(iatom)

    return iatoms

setattr(dump, "get_coord", get_coord)
setattr(dump, "get_iatom_from_type", get_iatom_from_type)
setattr(dump, "get_iatom_from_molid", get_iatom_from_molid)

if __name__ == '__main__':

    # parser arguments
    parser = Parser()
    args = parser.args

    # initialize parameters
    global para
    para = Parameters(args)

    # read masses
    global masses
    masses = [0.0]
    masses.extend(np.vectorize(float)(args.masses.split(',')))

    # mapping from special_atomid to its true value
    global special_atomids_values
    special_atomids_values = {}

    # read the input trj and select the frame we want
    global d, natoms, dof, timestep_index
    d = dump(args.input_lammpstrj)
    vx_indx = d.names['vx']; vy_indx = d.names['vy']; vz_indx = d.names['vz']
    timestep = int(args.timestep)
    if timestep == -1:
        timestep = d.time()[-1]
        timestep_index = -1
    if timestep in d.time():
        timestep_index = d.time().index(timestep)
    else:
        raise RuntimeError('cannot find timestep ' + str(timestep))
    d.tselect.one(timestep)
    print("selected %d in input trj"%timestep)
    natoms = d.snaps[timestep_index].natoms
    dof = 3 * natoms - 3

    global pbc
    snap = d.snaps[timestep_index]
    pbc = Pbc(snap)

    # parse hb_measure
    global hb_measure
    hb_measure = HB_Measure(args.hb_measure)

    # read the mode file to know which momentums are going to be kicked
    mode = Mode(args.mode)

    # parse coord_system
    coord_system = CoordSystem(
        args.coord_system, do_unwrapping = not args.donot_unwrapping)

    # construct coord system (this must be done after pbc is defined)
    coord_system.construct()
    if not args.quiet: coord_system.show()

    # decompose Mode in coord_systetm onto xyz
    mode.compute(coord_system, float(args.dT), kick_type = args.kick_type)
    if not args.quiet: mode.show()

    # apply the kick
    for iatom, movement in zip(mode.true_ids, mode.true_movements):
        snap.atoms[iatom - 1, [vx_indx, vy_indx, vz_indx]] += movement

    # write out resulting lammpstrj
    d.write(args.output_lammpstrj)
