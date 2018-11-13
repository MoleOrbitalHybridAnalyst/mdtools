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
import multiprocessing as mp
from pyemma.thermo import dtram,wham

class MyParser:

    def __init__(self):
        self.parser = argparse.ArgumentParser(\
                formatter_class = argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument('-m', '--metafile',
                required = True,
                help = 'input metadatafile')
        self.parser.add_argument('-n', '--cvname',
                default = 'cv',
                help = 'colvar name in time series')
        self.parser.add_argument('-T', '--temperature',
                required = True,
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
                help = 'job type, full or bootstrap or blockN')
        self.parser.add_argument('-l', '--lagtime',
                required = True,
                help = 'comma seperated lag times')
        self.parser.add_argument('-i', '--input',
                required = True,
                help = 'a json file controlling histo and bootstrap')
        self.parser.add_argument('--nprocs',
                default = 1,
                help = 'number of processes to be used')
        self.args = self.parser.parse_args()

class MetaData:

    def __init__(self, fname, kBT):
        # read metafile 
        # construct window_centers K filenames
        # window_centers[i] = center of i-th window
        # K[i] = force const of i-th in kT
        # filenames[i] = filename of i-th window
        self.filenames = []
        self.window_centers = []
        self.K = []
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
        self.nwindows = len(self.window_centers)
        self.nthermo = self.nwindows

    def bias_energy(self, thermo_index, x):
        return self.K[thermo_index]*(x-self.window_centers[thermo_index])**2

class MyDtram:

    def __init__(self, args):
        self.args = args
        self.metadata = MetaData(self.args.metafile, 
                      float(args.temperature) * float(args.kb))

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
        self.bias = []
        self.histo_centers = np.arange( \
                self.histo_min + self.binsize / 2,
                self.histo_max + self.binsize, self.binsize)
        for thermo_index, fn_cv in enumerate(self.metadata.filenames):
            self.dtrajs.append(self.make_dtraj(fn_cv))
            self.ttrajs.append(
                    np.ones(len(self.dtrajs[-1]), dtype = int) * thermo_index)
            self.bias.append(
                self.metadata.bias_energy(thermo_index, self.histo_centers))
        self.bias = np.array(self.bias)

        self.nprocs = int(self.args.nprocs)
        self.lags = np.vectorize(int)(self.args.lagtime.split(','))

    def make_dtraj(self, fname):
        df = read_colvar(fname)

        # get time info
        min_time = df.time.min()
        max_time = df.time.max()
        timestep = df.time[1] - df.time[0]
        length = int(round( (max_time - min_time) / timestep )) + 1 
        # make time series
        ts = np.ones(length) * np.inf;
        for time, cv in zip(df.time, df[self.args.cvname]):
            indx = int(round( (time - min_time) / timestep ))
            ts[indx] = cv

        where_is_missing = np.where(ts == np.inf)[0]
        where_is_present = np.where(ts != np.inf)[0]
        print("from %s read %d frames"%(fname,len(where_is_present)))
        missing_times =  where_is_missing * timestep + min_time
        for time in missing_times:
            print("missing time", time, "in", fname, file = stderr)

        # guess the missing values
        if len(where_is_missing):
            ts[where_is_missing] = \
                np.interp( \
                    missing_times, 
                    np.arange(
                        min_time, 
                        max_time + timestep, timestep)[where_is_present],
                    ts[where_is_present])

        #dtraj = np.digitize()
        #do not need digitize if binsizes are uniform
        dtraj = np.vectorize(int)( (ts - self.histo_min) / self.binsize)
        indx_max = int((self.histo_max - self.histo_min) / self.binsize)
        dtraj[dtraj < 0] = 0
        dtraj[dtraj > indx_max] = indx_max

        return dtraj

    def do_task(self, itask):
        if self.args.job_type == 'full':
            nlags = len(self.lags)
            if itask >= nlags: return None
            ilag = itask % nlags
            return [
                self.lags[ilag],
                float(self.args.temperature) * float(self.args.kb) * 
                dtram(
                    self.ttrajs, self.dtrajs, self.bias, 
                    self.lags[ilag], init = 'wham', 
                    init_maxiter = self.maxiter, init_maxerr = self.maxerr,
                    maxiter = self.maxiter, maxerr = self.maxerr).f_full_state
                ]
        elif self.args.job_type == 'bootstrap':
            pass
        elif re.match('block[0-9]+$', self.args.job_type):
            if itask >= len(self.lags) * self.nblocks: return None
            iblock = itask % self.nblocks
            ilag = int(itask / self.nblocks)
            blocksizes = [int(round(len(_)/self.nblocks)) for _ in self.ttrajs]
            if iblock != self.nblocks - 1:
                ttrajs = \
                  [t[iblock*s:(iblock+1)*s] for t,s in zip(self.ttrajs,blocksizes)]
                dtrajs = \
                  [t[iblock*s:(iblock+1)*s] for t,s in zip(self.dtrajs,blocksizes)]
            else:
                ttrajs = [t[iblock*s:] for t,s in zip(self.ttrajs,blocksizes)]
                dtrajs = [t[iblock*s:] for t,s in zip(self.dtrajs,blocksizes)]
            return [
                self.lags[ilag],
                float(self.args.temperature) * float(self.args.kb) * 
                dtram(
                    ttrajs, dtrajs, self.bias, 
                    self.lags[ilag], init = 'wham', 
                    init_maxiter = self.maxiter, init_maxerr = self.maxerr,
                    maxiter = self.maxiter, maxerr = self.maxerr).f_full_state
                ]
        else:
            raise Exception('unknown job type')

    def run(self, maxiter = 100000, maxerr = 1e-12):
        self.maxiter = maxiter
        self.maxerr = maxerr

        if self.args.job_type == 'full':
            ntasks = len(self.lags)
        elif self.args.job_type == 'bootstrap':
            ntasks = len(self.lags) * self.boot_nres
        elif re.match('block[0-9]+$', self.args.job_type):
            self.nblocks = int(self.args.job_type[5:])
            ntasks = len(self.lags) * self.nblocks
        else:
            raise Exception('unknown job type')

        with mp.Pool(int(self.args.nprocs)) as p:
            self.f_full_states = p.map(wrap_do_task, range(ntasks))

#        if self.args.job_type == 'full':
#            self.f_full_states = {
#                int(lag): 
#                    dtram( self.ttrajs, self.dtrajs, 
#                           self.bias, int(lag),
#                           init = 'wham',
#                           maxiter = self.maxiter, maxerr = self.maxerr
#                         ).f_full_state for lag in lags}
#        elif self.args.job_type == 'bootstrap':
#            pass
#        elif re.match('block[0-9]+$', self.args.job_type):
#            self.nblocks = int(self.args.job_type[5:])
#        else:
#            raise Exception('unknown job type')

    def save(self):
        prefix = self.args.output

        if self.args.job_type == 'full':
            for f in self.f_full_states:
                np.savetxt(
                    prefix + "_lag" + str(f[0]) + ".dat",
                    np.transpose([self.histo_centers, f[1]])
                )
        elif self.args.job_type == 'bootstrap':
            pass
        elif re.match('block[0-9]+$', self.args.job_type):
            iblock = -1
            for f in self.f_full_states:
                iblock = (iblock + 1) % self.nblocks
                np.savetxt(
                    prefix + "_lag" + str(f[0]) + 
                    "_block" + str(iblock) + ".dat",
                    np.transpose([self.histo_centers, f[1]])
                )
        else:
            raise Exception('unknown job type')

# ugly solution for using mp.Pool
def wrap_do_task(itask):
    return mydtram.do_task(itask)

if __name__=="__main__":

    parser = MyParser()
    mydtram = MyDtram(parser.args)
    mydtram.run(maxiter=500)
