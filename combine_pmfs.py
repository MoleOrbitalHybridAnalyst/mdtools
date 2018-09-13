# python3 combine_pmfs.py pmf0.dat pmf1.dat ...
from __future__ import print_function
import subprocess as sp
import numpy as np
from sys import argv
import re

nblocks = len(argv) - 1
cvs = []
pmfs = []
ndatas = []
for ifile,filename in enumerate(argv[1:]):
    pmf = []
    cv = []
    with open(filename, "r") as fp:
        for line in fp:
            if re.match("#", line)!=None:
                continue
            splits = line.split()
            cv.append(float(splits[0]))
            pmf.append(float(splits[1]))
    cvs.append(np.array(cv))
    pmfs.append(np.array(pmf))
    ndatas.append(len(pmf))
ndata = int(np.mean(ndatas))
cv_min = min([min(_) for _ in cvs])
cv_max = max([max(_) for _ in cvs])
cv_common = np.linspace(cv_min, cv_max, num=ndata, endpoint=True)
pmfs_interp = [np.interp(cv_common, _, __) for _,__ in zip(cvs, pmfs)]
averages = np.mean(pmfs_interp, axis = 0)
stddevs = np.std(pmfs_interp, axis = 0, ddof = 1)
for (cv, average, stddev) in zip(cv_common, averages, stddevs):
    print(cv,average,stddev)
