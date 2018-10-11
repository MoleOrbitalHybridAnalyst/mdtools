# python3 combine_pmfs.py pmf0.dat pmf1.dat ...
# assume files are aligned
# will generate pmf_averaged.dat and errors.dat
from __future__ import print_function
import subprocess as sp
import numpy as np
from sys import argv
import re

fn_output = "pmf_averaged.dat"
fn_errors = "errors.dat"

nblocks = len(argv) - 1
cvs = []
pmfs = []
ndatas = []
comments = []
for ifile,filename in enumerate(argv[1:]):
    print("reading file", filename)
    pmf = []
    cv = []
    if ifile == 0: line_count = 0
    with open(filename, "r") as fp:
        for line in fp:
            if re.match("#", line) != None or re.match("\s+$", line) != None:
                if ifile == 0:
                    comments.append([line_count, line])
                continue
            splits = line.split()
            cv.append([float(splits[0]), float(splits[1])])
            pmf.append(float(splits[2]))
            if ifile == 0: line_count += 1
    cvs.append(np.array(cv))
    pmfs.append(np.array(pmf))
    ndatas.append(len(pmf))

averages = np.mean(pmfs, axis = 0)
stddevs = np.std(pmfs, axis = 0, ddof = 1)

with open(fn_output, 'w') as fp_output:
    indx_comment = 0
    for indx, (cv, average) in enumerate(zip(cvs[0], averages)):
        if indx == comments[indx_comment][0]:
            print(comments[indx_comment][1], end = '', file = fp_output)
            indx_comment += 1
        print(cv[0], cv[1], average, file = fp_output)
print(fn_output, "generated")

with open(fn_errors, 'w') as fp_errors:
    indx_comment = 0
    for indx, (cv, stddev) in enumerate(zip(cvs[0], stddevs)):
        if indx == comments[indx_comment][0]:
            print(comments[indx_comment][1], end = '', file = fp_errors)
            indx_comment += 1
        print(cv[0], cv[1], stddev, file = fp_errors)
print(fn_errors, "generated")
