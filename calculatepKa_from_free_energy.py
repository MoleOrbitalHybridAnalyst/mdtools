#! /home/chhli/miniconda3/bin/python3
# read a free energy file
# calculate pKa and give the shifted potential of mean force
import argparse
import numpy as np
import re

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('free_data', help='input free energy file')
    parser.add_argument('--pmf_data', help='output pmf file')
    parser.add_argument('--ts_range', 
            help='find transition state in this range', default="1.7,2.7")
    parser.add_argument('--r_bulk', 
            help='consider to be bulk after this distance', default=8.0)
    parser.add_argument('--col', 
            help='column of the free energy in free.dat', default=1)
    parser.add_argument('--T', 
            help='temperature in Kelvin', default=310)
    parser.add_argument('--kb', 
            help='boltzmann constant in [free energy unit]/K', 
            default=0.0019872041)
    parser.add_argument('-v','--verbose', 
            help='print info to screen', 
            action="store_true")
    args = parser.parse_args()

    ts_ranges = args.ts_range.split(',')
    ts_ranges = np.vectorize(float)(np.array(ts_ranges))
    r_bulk = float(args.r_bulk)
    T = float(args.T)
    kb = float(args.kb)
    beta = 1/kb/T
    if r_bulk < ts_ranges[-1]:
        print("WARNING: r_bulk is smaller than the upper bound of ts_range")
    args.col = int(args.col)
    if args.verbose:
        print("kbT =",kb*T, "in energy unit")
        print("finding TS within %f and %f"%(ts_ranges[0],ts_ranges[1]))
        print("region beyond",r_bulk,"is regarded as bulk")

    free = []
    with open(args.free_data,"r") as fp:
        for line in fp:
            if re.match("#", line) != None:
                continue
            splits = line.split()
            free.append(
                    np.array([float(splits[0]),float(splits[args.col])]))
    free = np.array(free)
    if np.min(free[:,0])<0:
        print("ERROR: smallest r is",np.min(free[:,0]),"< 0")
        exit()
    if np.max(free[:,0])<args.r_bulk:
        print("ERROR: maximum r is",np.max(free[:,0]),"> r_bulk")
        exit()

    mask = free[:,0] > r_bulk
    ave = np.mean(free[mask,1]+kb*T*np.log(4*np.pi*free[mask,0]**2))
    free[:,1] -= ave
    if args.pmf_data != None:
        with open(args.pmf_data, "w") as fp:
            [print(x,y+kb*T*np.log(4*np.pi*x**2),file=fp) for (x,y) in free]
    mask = [x < ts_ranges[1] and x > ts_ranges[0] for x in free[:,0]]
    mask = np.array(mask)
    r_ts = free[mask,0][np.argmax(free[mask,1])]
    if args.verbose:
        print("transition state is found at",r_ts)
    
    mask = free[:,0] <= r_ts
    reduced = free[mask]
    integral = np.dot(free[1:,0]-free[:-1,0], 
            (np.exp(-beta*free[1:,1])+np.exp(-beta*free[:-1,1]))/2.0)
    print(np.log10(integral / 1660.0))
