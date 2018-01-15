# python3 dtram.py metafile temperature dilute_by
# metafile: 
# (each window each line; multi trjs permitted for each window)
# (use the force_const typically in used in plumed.dat)
#   file1 file2 ... center force_const
# dilute_by:
# (skipping every this values to reduce comp. cost)
from pyemma.thermo import dtram,wham
import numpy as np
import pandas as pd
from sys import argv

kB = 0.0019872041 # kB in kcal/mol

cvname = 'sol'
hist_start=0.2
hist_end=5.2
hist_resolution=0.025  # 200 bins
hist_centers=np.arange(hist_start,hist_end+hist_resolution,hist_resolution)

lags=[4,8]

nstates=len(hist_centers)
#timestep=0.5

# translate cv value into state index
def value2index(x):
    if((x - hist_start) / hist_resolution + 0.5 < 0): #if < min of bins
        conf_index=0
    else:
        conf_index=int( (x - hist_start) / hist_resolution + 0.5)
    if(conf_index >= nstates):
        conf_index=nstates-1
    return conf_index


def read_meta():
    # read metafile 
    # construct window_centers K filenames
    # window_centers[i] = center of i-th window
    # K[i] = force const of i-th in kT
    # filenames[i][j] = filename of j-th file of i-th window
    with open(argv[1], "r") as fp_meta:
        for line in fp_meta:
            splits = line.split()
            nfile = len(splits) - 2
            filenames.append(splits[:nfile])
            window_centers.append(float(splits[nfile]))
            # plumed use 0.5 * K * delta**2
            # assume K is in kcal/mol
            K.append(0.5 * float(splits[nfile + 1]) \
                    / float(argv[2]) / kB)

def bias_energy(thermo_index,y):
    return K[thermo_index]*(y-window_centers[thermo_index])**2

def make_traj(thermo_index,fn_cv):
    dtraj=[]
    ttraj=[]
    df = pd.read_csv(fn_cv)
    df.drop_duplicates(subset = 'time')
    full_len = df['time'].size
    min_time = df['time'].min()
    max_time = df['time'].max()
    chunk_size = (max_time - min_time) / full_len * float(argv[3])
    diluted_cvs = np.interp(np.arange(min_time, max_time, chunk_size), \
            df['time'].values, df[cvname].values)
    for cv in diluted_cvs:
        ttraj.append(thermo_index)
        dtraj.append(value2index(cv))
    print("from %s read %d frames"%(fn_cv,len(diluted_cvs)))
    return (ttraj,dtraj)

if __name__=="__main__":
    window_centers=[]; K=[]; filenames=[]
    bias=[]; dtrajs=[]; ttrajs=[]
    read_meta()
    thermo_index=0
    for thermo_index,cv in enumerate(window_centers):
        bias_row=[]
        for center in hist_centers:
            bias_row.append(bias_energy(thermo_index,center))
        bias.append(bias_row)
        for fn_cv in filenames[thermo_index]:
            (ttraj,dtraj)=make_traj(thermo_index,fn_cv)
            ttrajs.append(ttraj)
            dtrajs.append(dtraj)
    for thermo_index,cv in enumerate(window_centers):
        print("window %.3f  bias %f kT"%(cv,K[thermo_index]))
    bias=np.array(bias)
    dtram_objs=dtram(ttrajs,dtrajs,bias,lags,init='wham',maxiter=500000)
    for dtram_obj in dtram_objs:
    #    wham_obj=wham(ttrajs,dtrajs,bias)
        print();
        print(dtram_obj.f_full_state)
    #    print(wham_obj.f_full_state)
