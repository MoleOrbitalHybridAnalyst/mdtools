# reweight time series accordining to given reweighting factor by generating fake time series

from read_colvar import *
import argparse
from sys import stderr

def parse():
    parser = argparse.ArgumentParser(\
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_colvar', help = 'input colvar file')
    parser.add_argument('output_colvar', help = 'output colvar file')
    parser.add_argument('--histo', required = True, 
            help = 'histo_min,histo_max,nbins')
    parser.add_argument('--cv_name', default = 'cv', 
            help = 'cv name in  colvar file')
    parser.add_argument('--rf_name', default = 'rf', 
            help = 'reweighting factor name in colvar file ' +
            'or a file containing it')
    return parser.parse_args()

def read_rf():
    try:
        return df[args.rf_name]
    except:
        print('cannot find', args.rf_name, 'in', args.input_colvar, file = stderr)
        print('Try to read file named', args.rf_name, file = stderr)
    try:
        return np.loadtxt(args.rf_name)
    except Exception as e:
        print(e, file = stderr)
        exit(1)

if __name__ == '__main__':

    args = parse()

    # read in COLVAR
    try:
        df = read_colvar(args.input_colvar)
    except Exception as e:
        print(e, file = stderr)
        exit(1)

    # get cv from COLVAR
    try:
        cv_ts = df[args.cv_name]
    except:
        print('cannot find', args.cv_name, 'in', args.input_colvar, file = stderr)
        exit(1)

    # get reweighting factor
    rf = read_rf()

    # parse histo info
    histo_info = args.histo.split(',')
    if len(histo_info) != 3:
        print('incomplete histo info')
        exit(1)
    try:
        histo_min = float(histo_info[0])
        histo_max = float(histo_info[1])
        nbins = int(histo_info[2])
    except Exception as e:
        print('while parsing histo info, ', end = '', file = stderr)
        print(e, file = stderr)
        exit(1)

    # check things
    if histo_min >= histo_max:
        print('histo_min should be smaller than histo_max')
        exit(1)
    if len(cv_ts) != len(rf):
        print('length of cv time series (%d) !='%len(cv_ts),
              'length of reweighting factor (%d)'%len(rf))
        exit(1)

    # do histogram using rf
    counts = np.zeros(nbins)
    npoint = 0
    for cv, w in zip(cv_ts, rf):
        ibin = int((cv - histo_min) / (histo_max - histo_min) * nbins)
        if ibin >= 0 and ibin < nbins:
            counts[ibin] += w
            npoint += 1
    if sum(counts) == 0:
        print('no points lie in histo', file = stderr)
        exit(1)
    counts /= sum(counts) 
    counts = np.vectorize(round)(counts * npoint)
    counts = np.vectorize(int)(counts)
    if npoint != sum(counts):
        diff = npoint - sum(counts)
        print('WARNING: %d point(s) missing,'%diff,
              'estimated error = %f%%'%(100 * diff / npoint),
                file = stderr)
        maxpos = np.argmax(counts)
        counts[maxpos] += diff
        if counts[maxpos] <= 0 or np.argmax(counts) != maxpos:
            print('too much difference between two ensembles', file = stderr)
            exit(1)
