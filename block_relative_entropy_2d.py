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
    parser.add_argument('-n', '--cv_names', default = 'cv,sol', \
            help = 'input colvar names')
    parser.add_argument('-c', '--cv_columns', default = '1,2', \
            help = 'input colvar columns')
    parser.add_argument('-b', '--nblocks', default = 10, \
            help = 'number of blocks the input will be split into')
    parser.add_argument('-w','--widths', \
            help='widths of bins', required = True)
    parser.add_argument('-r','--ratio', \
            help='extend the histo range by this ratio of true range', \
            default = 0.1)

    return parser.parse_args()

def make_histo2d(array, binsizes, ratio, min_maxs = None):
    mask = np.sum(np.isfinite(array), axis = 1) == 2
    array = array[mask]
    if min_maxs is None:
        min_x = min(array[:, 0]); max_x = max(array[:, 0])
        extend = (max_x - min_x) * ratio / 2.0
        min_x -= extend; max_x += extend
        min_y = min(array[:, 1]); max_y = max(array[:, 1])
        extend = (max_y - min_y) * ratio / 2.0
        min_y -= extend; max_y += extend
        min_maxs = [[min_x, max_x], [min_y, max_y]]
    min_x = min_maxs[0][0]; max_x = min_maxs[0][1]
    min_y = min_maxs[1][0]; max_y = min_maxs[1][1]
    if debug:
        print('min_x =', min_x, 'max_x =', max_x)
        print('min_y =', min_y, 'max_y =', max_y)
    nbins_x = int(round((max_x - min_x) / binsizes[0]))
    nbins_y = int(round((max_y - min_y) / binsizes[1]))
    histo = np.zeros(nbins_x * nbins_y)
    for a in array:
        ix = min(int((a[0]-min_x)/binsizes[0]), nbins_x - 1)
        iy = min(int((a[1]-min_y)/binsizes[1]), nbins_y - 1)
        histo[ix * nbins_y + iy] += 1
    total = sum(histo)
    return [histo / total, min_maxs]

def relative_entropy(dist1, dist2, epsilon):
    """
        reference_dist, model_dist
    """
    assert len(dist1) == len(dist2)
    result = 0.0
    for x,y in zip(dist1, dist2):
        if x < 1e-100:
            continue
        if y < 1e-100:
            result += x * np.log(x/(y+epsilon))
            continue
        result += x * np.log(x/y)
    return result

if __name__ == "__main__":

    debug = 0

    args = parse()

    if args.file_format == 'colvar':
        print("# info: read " + args.input + " as colvar")
        if args.cv_names is None:
            print('cv_names is required', file = stderr); exit()
        cv_names = args.cv_names.split(',')
        if len(cv_names) != 2:
            print('two cv names expected', file = stderr); exit()
        print("# info: read fields %s %s"%(cv_names[0], cv_names[1]))
        from read_colvar import read_colvar
        time_series = \
            read_colvar(args.input)[cv_names].values
    elif args.file_format == 'data':
        print("# info: read " + args.input + " as data")
        try:
            cols = np.vectorize(int)(args.cv_columns.split(','))
        except:
            print('integer cv_columns is required', file = stderr); exit()
        if len(cols) != 2:
            print('two columns expected', file = stderr); exit()
        time_series = np.loadtxt(args.input)[:, cols]
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

    widths = np.vectorize(float)(args.widths.split(','))
    for width in widths:
        if width <= 0.0:
            print('width should be positive', file = stderr); exit()

    full_histo, full_min_max = make_histo2d(time_series, widths, ratio)
    if debug:
        print("full_histo")
        print(full_histo)

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
        histo, junk =  \
            make_histo2d(time_series[lo:hi], widths, ratio, min_maxs = full_min_max)
        if debug:
            print("histo",iblock)
            print(histo)
        print(iblock, relative_entropy(full_histo, histo, epsilon))
