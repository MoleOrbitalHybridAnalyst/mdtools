# assume mol type x y z vx vy vz are in input_lammpstrj

# mode file format:
#  id0 dx0 dy0 dz0
#  id1 dx1 dy1 dz1
#  id2 dx2 dy2 dz2

# note that using special atom ids are assuming there is only one proton
# behaviour of using special atom ids with multiple protons is undefined

# special atom ids:
# OTH: atom whose type == type_hyd_o
# OW1: water oxygen which has strongest HB with OTH
# OW2: water oxygen which has second strongest HB with OTH
# OW3: water oxygen which has third strongest HB with OTH
# HW1*: water hydrogen which is bonded to OW1
# HW2*: water hydrogen which is bonded to OW2
# HW3*: water hydrogen which is bonded to OW3

from lammpstrj import *
from argparse import *

class Parameters:
    """ hard coded parameters used in this code """

    def __init__(self, args):
        self.args = args

        self.oxygen_search_cutoff = 6.0

        self.special_atomids = \
              [ 'OTH' , 'HTH1', 'HTH2', 'HTH3', \
                'OW1' , 'OW2' , 'OW3' ,         \
                'HW11', 'HW21', 'HW31',         \
                'HW12', 'HW22', 'HW32'  ]
        try: self.type_water_o = int(args.type_water_o)
        except TypeError: pass
        try: self.type_water_o = int(args.type_water_o)
        except TypeError: pass
        try: self.type_hyd_o = int(args.type_hyd_o)
        except TypeError: pass
        try: self.type_hyd_h = int(args.type_hyd_h)
        except TypeError: pass

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
        parser.add_argument('--hb_measure', '-h',
                help = 'descriptor of hydrogen bonds',
                default = None)
        parser.add_argument('--do_unwrapping', '-u',
                help = 'do unwrapping when defining coordinate system'
                action = 'store_true' )
        parser.add_argument('--type_water_o', help = 'type of water oxygen')
        parser.add_argument('--type_water_h', help = 'type of water hydrogen')
        parser.add_argument('--type_hyd_o', help = 'type of hydronium oxygen')
        parser.add_argument('--type_hyd_h', help = 'type of hydronium hydrogen')
        self.args = parser.parse_args()

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
        global pbc
        measure_para = hb_measure_string.split(',')
        if len(hb_measure_string) == 0:
            return None
        self.measure_type = measure_para[0]
        # parse measure_type
        if measure_type == 'delta':
            if len(measure_para) != 1:
                raise RuntimeError('no parameters needed for delta')
        elif measure_type == 'roo_delta':
            if len(measure_para) != 3:
                raise RuntimeError('roo and theta are needed for roo_delta')
            self.roo = float(measure_para[1])
            self.theta = float(measure_para[2])
        elif measure_type == 'roo_angle':
            raise RuntimeError('not implemented yet')
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
        else:
            raise RuntimeError('cannot understand hb_measure', measure_type)

    def translate_atomid(atomid, para, lmptrj, special_atomids_values):
        if type(atomid) is int: return atomid
        elif type(atomid) is str:
            if atomid in para.special_atomids:
                pass
            else:
                raise RuntimeError('cannot understand atomid ' + atomid)
    
            if atomid[0:3] == 'OTH':
        else:
            raise RuntimeError('int or str expected in translate_atomid')


class CoordSystem:
    """ class of coordinate system """

    def __init__(self, coord_system_string, \
            do_unwrapping = True, hb_measure = None):
        self.hb_measure = hb_measure
        self.do_unwrapping = do_unwrapping
        self.ex = np.array([1.0, 0.0, 0.0])
        self.ey = np.array([0.0, 1.0, 0.0])
        self.ez = np.array([0.0, 0.0, 1.0])
        #--- read atom ids for coodinate system setup --#
        self.ids = coord_system_string.split(',')
        if len(self.ids) != 3:
            raise RuntimeError('3 atom ids expected for coordinate system')

    def construct():
        #TODO if HTO HTH OW1 OW2 OW3 are used, hb_measure must be set
        # if x,y,z, then pass
        # PBC !!!
        # if x,y,z are mixed with id0, id1, id2, raise Error
        # if id0, id1, id2, then ex is along id0->id2; ez is id0->id2 x id0->id1
        # ey is ez x ex
        # if ez === \vec{0}, ez = (0, 0, 1.0) - ((0, 0, 1).ex) ex
        global hb_measure
        global para
        for atomid in ids:
            if atomid in para.special_atomids:
                if 
                # TODO translate special atomid

class Mode:
    """ class of motion mode """

    def __init__(self, para):
        data = np.genfromtxt(para.mode, dtype = str)
        self.ids = [line[0] for line in data]
        self.movements = [np.vectorize(float)(line[1:4]) for line in data]

def get_coord(lmptrj, iatom):
    return lmptrj.atom(iatom, 'x', 'y', 'z').flatten()

setattr(dump, "get_coord", get_coord)

if __name__ == '__main__':

    #--- parser arguments ---#
    parser = Parser()
    args = parser.args

    #--- initialize parameters ---#
    global para
    para = Parameters(args)

    # mapping from special_atomid to its true value
    global special_atomids_values
    special_atomids_values = {}

    #--- parse HB_Measure ---#
    global hb_measure
    hb_measure = HB_Measure(args.hb_measure)

    #--- read the mode file to know which momentums are going to be kicked ---#
    mode = Mode(para)

    # read the input trj
    d = dump(args.input_lammpstrj)

    # select the frame we want
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

    snap = d.snaps[timestep_index]

    global pbc
    pbc = Pbc(snap)
