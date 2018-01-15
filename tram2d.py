# python3 tram2d.py metafile temperature dilute_by
# metafile: 
# (each window each line; multi trjs permitted for each window)
# (use the force_const typically in used in plumed.dat)
# (markov model will be computed on window in the first line)
#   file1 file2 ... center force_const
# dilute_by:
# (skipping every this values to reduce comp. cost)
from pyemma.thermo import tram,wham
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eig
import pandas as pd
from sys import argv

kB = 0.0019872041 # kB in kcal/mol

# cv 1
cvname1 = 'rmsd'
hist_start1=0.042
hist_end1=0.1
hist_resolution1=0.0116  # 5 bins
hist_centers1=np.arange(hist_start1,hist_end1+hist_resolution1,hist_resolution1)
# cv 2
cvname2 = 'sol'
hist_start2=0.2
hist_end2=5.2
hist_resolution2=0.1    #  50 bins
hist_centers2=np.arange(hist_start2,hist_end2+hist_resolution2,hist_resolution2)

lags=[1,2,4,8]

nstates1=len(hist_centers1)
nstates2=len(hist_centers2)
#timestep=0.5

# translate cv1 value into state index
def value2index1(x):
    if((x - hist_start1) / hist_resolution1 + 0.5 < 0): #if < min of bins
        conf_index=0
    else:
        conf_index=int( (x - hist_start1) / hist_resolution1 + 0.5)
    if(conf_index >= nstates1):
        conf_index=nstates1-1
    return conf_index

# translate cv2 value into state index
def value2index2(x):
    if((x - hist_start2) / hist_resolution2 + 0.5 < 0): #if < min of bins
        conf_index=0
    else:
        conf_index=int( (x - hist_start2) / hist_resolution2 + 0.5)
    if(conf_index >= nstates2):
        conf_index=nstates2-1
    return conf_index

# translate cv value pair into state index
def value2index(x,y):
    i1 = value2index1(x)
    i2 = value2index2(y)
    return i1 * nstates2 + i2

def read_meta():
    # read metafile 
    # construct window_centers K filenames
    # window_centers1[i] = center of i-th window of cv1
    # K1[i] = force const of i-th window in kT
    # filenames[i][j] = filename of j-th file of i-th window
    with open(argv[1], "r") as fp_meta:
        for line in fp_meta:
            splits = line.split()
            nfile = len(splits) - 4
            filenames.append(splits[:nfile])
            window_centers1.append(float(splits[nfile]))
            window_centers2.append(float(splits[nfile + 2]))
            # plumed use 0.5 * K * delta**2
            # assume K is in kcal/mol
            K1.append(0.5 * float(splits[nfile + 1]) \
                    / float(argv[2]) / kB)
            K2.append(0.5 * float(splits[nfile + 3]) \
                    / float(argv[2]) / kB)

def bias_energy(thermo_index,x,y):
    return K1[thermo_index]*(x-window_centers1[thermo_index])**2 \
            + K2[thermo_index]*(y-window_centers2[thermo_index])**2 \

def make_traj(thermo_index,fn_cv):
    dtraj=[]; ttraj=[]; btraj=[]
    df = pd.read_csv(fn_cv)
    df.drop_duplicates(subset = 'time')
    full_len = df['time'].size
    min_time = df['time'].min()
    max_time = df['time'].max()
    chunk_size = (max_time - min_time) / full_len * float(argv[3])
    diluted_cvs1 = np.interp(np.arange(min_time, max_time, chunk_size), \
            df['time'].values, df[cvname1].values)
    diluted_cvs2 = np.interp(np.arange(min_time, max_time, chunk_size), \
            df['time'].values, df[cvname2].values)
    for (cv1,cv2) in zip(diluted_cvs1,diluted_cvs2):
        ttraj.append(thermo_index)
        dtraj.append(value2index(cv1,cv2))
        #btraj[i][j] is the bias energy of i-th frame evaluated in j-th thermostate
        #btraj will be btraj_elements
        btraj_element = [] 
        for thermo_index2 in range(len(window_centers1)):
            btraj_element.append(bias_energy(thermo_index2,cv1,cv2))
        btraj.append(np.array(btraj_element))
    print("from %s read %d frames"%(fn_cv,len(diluted_cvs1)))
    return (ttraj,dtraj,btraj)

if __name__=="__main__":
    window_centers1=[]; K1=[]; filenames=[]
    window_centers2=[]; K2=[]; filenames=[]
    btrajs=[]; dtrajs=[]; ttrajs=[]
    read_meta()
    for thermo_index,(cv1,cv2) in enumerate(zip(window_centers1,window_centers2)):
        for fn_cv in filenames[thermo_index]:
            (ttraj,dtraj,btraj)=make_traj(thermo_index,fn_cv)
            ttrajs.append(ttraj); dtrajs.append(dtraj); btrajs.append(btraj)
        print("%s: center = %.3f bias = %f kT; %s: center = %.3f bias = %f kT"\
                %(cvname1,cv1,K1[thermo_index],cvname2,cv2,K2[thermo_index]))
    #for thermo_index,(cv1,cv2) in enumerate(zip(window_centers1):
    #    print("window %.3f  bias %f kT"%(cv,K[thermo_index]))
    dtram_objs=tram(ttrajs,dtrajs,np.array(btrajs),lags,unbiased_state=0,init='wham',maxiter=500000,maxerr=1e-14)
    for dtram_obj in dtram_objs:
    #    wham_obj=wham(ttrajs,dtrajs,bias)
        print();
        print(dtram_obj.f_full_state)
    #    print(wham_obj.f_full_state)

    # past processing
    P = dtram_objs[-1].msm.P
    eigsys = eig(P)
    order = np.argsort(-abs(eigsys[0]))
    slowest = transpose(eigsys[1])[order[1]]
    with open("slowest.dat","w") as fp:
        for i1, cv1 in enumerate(hist_centers1):
            for i2, cv2 in enumerate(hist_centers2):
                print(cv1, cv2, slowest[i1*nstates+i2],file=fp)
            print("",file=fp)
