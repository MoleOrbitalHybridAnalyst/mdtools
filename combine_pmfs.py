# python3 combine_pmfs.py pmf0.dat pmf1.dat ...
import subprocess as sp
import numpy as np
from sys import argv
import re

nblocks = len(argv) - 1
cvs = []
pmfs = []
for ifile,filename in enumerate(argv[1:]):
    pmf = []
    with open(filename, "r") as fp:
        for line in fp:
            if re.match("#", line)!=None:
                continue
            splits = line.split()
            if ifile == 0:
                cvs.append(float(splits[0]))
            pmf.append(float(splits[1]))
    pmfs.append(np.array(pmf))
averages = np.mean(pmfs, axis = 0)
stddevs = np.std(pmfs, axis = 0, ddof = 1)
for (cv, average, stddev) in zip(cvs, averages, stddevs):
    print(cv,average,stddev)
