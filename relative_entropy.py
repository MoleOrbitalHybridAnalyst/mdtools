# calculate the relative entropy by giving two time series
import argparse
import re
import numpy as np

def make_histos(file1, file2):
    cv_trj1 = []; cv_trj2 = []
    with open(file1, "r") as fp:
        for line in fp:
            if re.match("#", line) != None:
                continue
            splits = line.split()
            cv_trj1.append(float(splits[c1]))
    with open(file2, "r") as fp:
        for line in fp:
            if re.match("#", line) != None:
                continue
            splits = line.split()
            cv_trj2.append(float(splits[c2]))
    cv_trj1 = np.array(cv_trj1)
    cv_trj2 = np.array(cv_trj2)
    min_ = min(cv_trj1); max_ = max(cv_trj1)
    min_ = min(min_,min(cv_trj2)); max_ = max(max_,max(cv_trj2))
    extend = (max_ - min_) * r / 2.0
    min_ -= extend; max_ += extend
    nbins = int((max_ - min_) / width) + 1
    histo1 = np.zeros(nbins)
    histo2 = np.zeros(nbins)
    for _ in cv_trj1:
        histo1[int((_-min_)/width)] += 1 
    for _ in cv_trj2:
        histo2[int((_-min_)/width)] += 1 
    histo1 /= len(cv_trj1)
    histo2 /= len(cv_trj2)
    return [len(cv_trj1), histo1, len(cv_trj2), histo2]

def relative_entropy(dist1, dist2, epsilon):
    result = 0.0
    for x,y in zip(dist1, dist2):
        if x < 1e-100:
            continue
        if y < 1e-100:
            result -= x * np.log(y+epsilon)
            continue
        result -= x * np.log(y)
    return result

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cv_file1', help='the cv file used as the target distribution')
    parser.add_argument('cv_file2', help='the cv file used as the model distribution')
    parser.add_argument('-c','--column',
            help='column numbers of cv in file1 and file2',default="1,1")
    parser.add_argument('-w','--width',
            help='width of bins in when doing histogram',default=0.025)
    parser.add_argument('-r','--ratio',
            help='extend the histo range by this ratio of true range',
            default=0.2)
    args = parser.parse_args()

    c1 = int(args.column.split(',')[0])
    c2 = int(args.column.split(',')[1])
    width = float(args.width)
    r = float(args.ratio)

    l1,h1,l2,h2 = make_histos(args.cv_file1, args.cv_file2)
    print(relative_entropy(h1, h2, 1.0/l2))

