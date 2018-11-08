# python3 dtram.py metafile temperature dilute_by
# metafile: 
# (each window each line; multi trjs permitted for each window)
# (use the force_const typically in used in plumed.dat)
#   file1 file2 ... center force_const
#
# bootstrap:
#   nresampling: number of resamplings to be done
#   chunksize:   ratio of size of chunk over lag time
# histo:
#   min:     histo lower boundary
#   max:     histo upper boundary
#   nbins:   number of bins
#   binsize: binsize
import re
import json
import argparse
from sys import stderr
from read_colvar import *
from pyemma.thermo import dtram,wham

class MyParser:

    def __init__(self):
        self.parser = argparse.ArgumentParser(\
                formatter_class = argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument('-m', '--metafile',
                help = 'input metadatafile')
        self.parser.add_argument('-n', '--cvname',
                default = 'cv',
                help = 'colvar name in time series')
        self.parser.add_argument('-T', '--temperature',
                help = 'temperature used in simulations')
        self.parser.add_argument('--dilute_by',
                help = 'dilute time series by this value to ' +
                       'save computational cost')
        self.parser.add_argument('--kb',
                default = 0.0019872041,
                help = 'Boltzmann constant')
        self.parser.add_argument('-o', '--output',
                default = 'pmf',
                help = 'output file prefix')
        self.parser.add_argument('-j', '--job_type',
                default = 'full',
                help = 
                'job type, full or bootstrap or blockN')
        self.parser.add_argument('-l', '--lagtime',
                help = 'comma seperated lag times')
        self.parser.add_argument('-i', '--input',
                help = 'a json file controlling histo and bootstrap')
        self.args = self.parser.parse_args()

class MetaData:

    def __init__(self, fname, kBT):
        # read metafile 
        # construct window_centers K filenames
        # window_centers[i] = center of i-th window
        # K[i] = force const of i-th in kT
        # filenames[i] = filename of i-th window
        with open(fname, "r") as fp_meta:
            for line in fp_meta:
                if re.match("\s*#", line): continue
                if re.match("\s*$", line): continue
                splits = line.split()
                self.filenames.append(splits[0])
                self.window_centers.append(float(splits[1]))
                # plumed use 0.5 * K * delta**2
                # covert K into kBT
                self.K.append(0.5 * float(splits[2]) / kBT)
        self.nwindows = len(window_centers)
        self.nthermo = self.nwindows

    def bias_energy(self, thermo_index, x):
        return self.K[thermo_index]*(x-self.window_centers[thermo_index])**2

class MyDtram:

    def __init__(self, args):
        self.args = args
        self.metadata = MetaData(self.args.metafile, 
                      float(args.T) * float(args.kb))

        fp = open(self.args.input)
        self.J = json.load(fp)
        fp.close()

        # parse histo
        if 'histo' in self.J:
            J_histo = self.J['histo']
            if 'min' in J_histo:
                self.histo_min = float(J_histo['min'])
            else:
                raise Exception('cannot find histo min setting')
            if 'max' in J_histo:
                self.histo_max = float(J_histo['max'])
            else:
                raise Exception('cannot find histo max setting')
            if 'nbins' in J_histo:
                nbins = int(J_histo['nbins'])
                self.binsize = (self.histo_max - self.histo_min) / nbins
            elif 'binsize' in J_histo:
                self.binsize = float(J_histo['binsize'])
            else:
                raise Exception('cannot find histo nbins or binsize setting')
        else:
            raise Exception('cannot find histo settings')

        # parse bootstrap
        if 'bootstrap' in self.J:
            J_boot = self.J['bootstrap']
            if 'chunksize' in J_boot:
                self.boot_chunksize = int(J_boot['chunksize'])
            else:
                raise Exception('cannot find bootstrap chunksize setting')
            if 'nresampling' in J_boot:
                self.boot_nres = int(J_boot['nresampling'])
            else:
                raise Exception('cannot find bootstrap nresampling setting')
        else:
            raise Exception('cannot find bootstrap settings')

        self.dtrajs = []
        self.ttrajs = []
        self.btrajs = []
        self.histo_centers = np.arange(self.binsize / 2, \
                     self.histo_max + self.binsize / 2, self.binsize)
        for thermo_index, fn_cv in enumerate(self.metadata.filenames):
            self.dtrajs.append(self.make_dtraj(fn_cv))
            self.ttrajs.append(
                    np.ones(len(self.dtrajs[-1]), dtype = int) * thermo_index)
            self.btrajs.append(
                self.metadata.bias_energy(thermo_index, self.histo_centers))

    def make_dtraj(self, fname):
        df = read_colvar(fname)

        # get time info
        min_time = df.time.min()
        max_time = df.time.max()
        timestep = df.time[1] - df.time[0]
        length = int(round( (max_time - min_time) / timestep )) + 1

        # make time series
        ts = np.ones(length) * np.inf;
        for time, cv in zip(df.time, df[self.args.cv]):
            indx = int(round( (time - min_time) / timestep ))
            ts[indx] = cv

        #dtraj = np.digitize()
        #do not need digitize if binsizes are uniform
        dtraj = np.vectorize(int)( (ts - self.histo_min) / self.binsize)
        indx_max = int((self.histo_max - self.histo_min) / self.binsize)
        dtraj[dtraj < 0] = 0
        dtraj[dtraj > indx_max] = indx_max
        print("from %s read %d frames"%(fname,len(dtraj)))

        # print warnings about missig times
        where_is_missing = np.where(ts == np.inf)[0]
        missing_times =  where_is_missing * timestep + min_time
        for time in missing_times:
            print("missing time", time, "in", fname, file = stderr)
        for i in where_is_missing:
            cv_guess = (ts[i-1] + ts[i+1]) / 2.0
            dtraj[i] = int( (cv_guess - self.histo_min) / self.binsize)

        return dtraj

if __name__=="__main__":

    parser = MyParser()
    dtram = MyDtram(parser.args)
