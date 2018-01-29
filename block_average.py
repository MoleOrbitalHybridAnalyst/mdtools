# do block average for one cv time series
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_colvar(fname):
    with open(fname, "r") as fp:
        firstline = fp.readline()
    if re.match("#!.+FIELDS", firstline) == None:
        raise NameError('No header found in ' + fname)
    names_ = firstline.split()[2:]
    return pd.read_csv(fname,delim_whitespace=True,names=names_, comment ='#' )

def block_average(values, nblock):
    stride = int(len(values)/nblock)
    indexes = np.arange(0, nblock * stride + 1, stride)
    slices = []
    for i in range(len(indexes)-1):
        slices.append(np.mean(values[indexes[i]:indexes[i+1]]))
    slices = np.array(slices)
    return [np.mean(slices),np.sqrt(np.var(slices,ddof=1)/len(slices))]

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cv_file', help='input colvar file')
    parser.add_argument('--cv_name', 
            help='cv name in colvar file', default="cv")
    parser.add_argument('-m', '--min_nblock',
            help='minimum number of blocks', default="3")
    args = parser.parse_args()

    df = read_colvar(args.cv_file)
    values = df[args.cv_name].values
    min_block = int(args.min_nblock)
    
    aves = []; stds = []; nblocks = []
    nblock = len(values)
    while nblock >= min_block:
        ave, std = block_average(values, nblock)
        aves.append(ave); stds.append(std); nblocks.append(nblock)
        nblock = int(nblock / 2)
    plt.plot(nblocks, stds)
    plt.show()
