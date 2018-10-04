"""
    read a colvar / data
    compute relative_entropy of each block compared to the whole file
"""
import re
import argparse
import numpy as np
from sys import stderr

def parse():

    parser = argparse.ArgumentParser( \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('input', help = 'input time series')
    parser.add_argument('-f', '--file_format', default = 'colvar', \
            help = 'input colvar file format (data or colvar)')
    parser.add_argument('-n', '--cv_name', default = 'cv', \
            help = 'input colvar name')
    parser.add_argument('-c', '--cv_column', default = 1, \
            help = 'input colvar column')
    parser.add_argument('-n', '--nblocks', default = 10, \
            help = 'number of blocks the input will be split into')
    parser.add_argument('-w','--width', \
            help='width of bins', required = True)
    parser.add_argument('-r','--ratio', \
            help='extend the histo range by this ratio of true range', \
            default = 0.1)

    return parser.parse_args()

def make_histo(array, binsize, ratio, min_max = None):
    if min_max is None:
        min_ = min(array); max_ = max(array)
        extend = (max_ - min_) * ratio / 2.0
        min_ -= extend; max_ += extend
        min_max = [min_, max_]
    min_ = min_max[0]; max_ = min_max[1]
    nbins = int((max_ - min_) / binsize)
    histo = np.zeros(nbins)
    for a in array:
        histo[int((a-min_)/width)] += 1 
    return [histo, min_max]

def relative_entropy(dist1, dist2, epsilon):
    """
        reference_dist, model_dist
    """
    result = 0.0
    for x,y in zip(dist1, dist2):
        if x < 1e-100:
            continue
        if y < 1e-100:
            result += x * np.log(x/(y+epsilon))
            continue
        result += x * np.log(x/y)
    return result

if __name__=="__main__":

    args = parse()

    if args.file_format == 'colvar':
        print("# info: read " + args.input + " as colvar")
        if args.cv_name is None:
            print('cv_name is required', file = stderr); exit()
        print("# info: read field " + args.cv_name)
        from read_colvar import read_colvar
        time_series = read_colvar(args.input)[args.cv_name].values
    elif args.file_format == 'data':
        print("# info: read " + args.input + " as data")
        try:
            col = int(args.cv_column)
        except:
            print('cv_column is required', file = stderr); exit()
        time_series = np.loadtxt(args.input)[:, col]
    else:
        print('unkown file_format', file = stderr); exit()

    ndata = len(time_series)
    print("# info: total number of data points read:", ndata)

    try:
        ratio = float(args.ratio)
    except:
        print('ratio should be a float', file = stderr); exit()
    if ratio <= 0.0:
        print('ratio should be positive', file = stderr); exit()

    width = float(args.width)
    if width <= 0.0:
        print('width should be positive', file = stderr); exit()

    full_histo, full_min_max = make_histo(time_series, width, ratio)

    try:
        nblocks = int(args.nblocks)
    except:
        print('nblocks should be an int', file = stderr); exit()
    if nblocks <= 0:
        print('nblocks should be positive', file = stderr); exit()

    block_size = int(ndata / nblocks)
    print("# info: number of data points in each block:", block_size)
    # minimum flutuation in prob. of block data
    epsilon = 1.0 / block_size

    for iblock in range(nblocks):
        lo = iblock * block_size
        hi = (iblock + 1) * block_size
        histo, =  \
        make_histo(time_series[lo, hi], width, ratio, min_max = full_min_max)
        print(iblock, relative_entropy(full_histo, histo, epsilon))
