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

    # check things
    if len(cv_ts) != len(rf):
        print('length of cv time series (%d) !='%len(cv_ts),
              'length of reweighting factor (%d)'%len(rf))
        exit(1)
