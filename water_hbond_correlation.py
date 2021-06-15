# python this.py input.xyz boxsize hb_corr.dat
import MDAnalysis as mda
from MDAnalysis.coordinates.XYZ import XYZReader
from sys import argv
from autocorrelate import *
from time import time

rcut = 3.5   # in angstrom
acut = 30    # in degree

# read xyz and gather some basic info
xyz = XYZReader(argv[1])
box = float(argv[2])
box = np.array([box]*3)
natoms = xyz.n_atoms
assert natoms % 3 == 0
nwat = natoms // 3
types = list()
with open(xyz.filename) as f:
    for i, l in enumerate(f):
        if i < 2: 
            continue
        if i == natoms + 2:
            break
        types.append(l.split()[0])
types = np.array(types)
print(f"{nwat} waters read")
oindexes = np.where(types == 'O')[0]
hindexes = np.where(types == 'H')[0]
N_O = len(oindexes)
N_H = len(hindexes)
assert N_O == nwat
assert N_H == 2 * nwat
rangeO = np.arange(N_O)
rangeH = np.arange(N_H)
rcut2 = rcut**2
cosacut = np.cos(acut*np.pi/180)


hbond_mat_trj = list()
nframes = xyz.n_frames
tenth_nframes = int(nframes/10)

for iframe, ts in enumerate(xyz):
    # N x M x 3
    doh = ts.positions[oindexes][:,None,:] - ts.positions[hindexes][None,...] 
    n = (doh + 0.5 * box) // box
    doh = doh - n * box
    doo = ts.positions[oindexes][:,None,:] - ts.positions[oindexes][None,...] 
    n = (doo + 0.5 * box) // box
    doo = doo - n * box

    # N x M
    doh2 = np.sum(doh**2, axis = 2)
#    doo2 = np.sum(doo**2, axis = 2)

    if iframe == 0:
        # build bond topo:

        # argsort among hydrogens
        # the bonded H should be the first two
        bondedHList = np.argsort(doh2, axis=1)[:,:2]

        check_h = np.unique(bondedHList.flatten())
        assert len(check_h) == N_H
        assert (np.unique(check_h) == rangeH).all()

    # the closest non-boned O to each H
    # argsort among oxygens 
    # need to determine if 0-th or the 1st is the bonded O
    nearestNBOList = np.argsort(doh2, axis=0)[:2] 
    # check if 0th oxygen is the bonded O
    check_h = bondedHList[nearestNBOList[0]]
    notBondedTo0 = \
        (check_h[:,0] != rangeH) & (check_h[:,1] != rangeH)
    # check if 1st oxygen is the bonded O
    check_h = bondedHList[nearestNBOList[1]]
    notBondedTo1 = \
        (check_h[:,0] != rangeH) & (check_h[:,1] != rangeH)
    # every H should be bonded to closest O or next closest
    assert ( np.array(notBondedTo0,dtype=int) \
           + np.array(notBondedTo1,dtype=int) == np.ones(N_H) ).all()
    nearestNBOList = nearestNBOList[np.array([notBondedTo0,notBondedTo1])]

#    distances = list()
#    indexes = list()
#    for i in range(natoms):
#        indx = np.argsort(d2ij[i])[:3]
#        d2 = d2ij[i][indx]
#        indexes.append(indx)
#        distances.append(np.sqrt(d2))

    hbond_mat = np.zeros([nwat, nwat])
    #       O_0
    #      /   \
    #     H_1  H_2
    #     .     .
    #    .       .
    #   O_1      O_2
    H1 = bondedHList[:,0]
    H2 = bondedHList[:,1]
    O1 = nearestNBOList[H1]
    O2 = nearestNBOList[H2]
    d2OO1 = np.sum(doo[rangeO, O1]**2, axis=1)
    d2OO2 = np.sum(doo[rangeO, O2]**2, axis=1)
    hbond_mat[rangeO, O1] += (d2OO1 < rcut2)
    hbond_mat[rangeO, O2] += (d2OO2 < rcut2)
    vOH1 = doh[rangeO, H1]
    vOH2 = doh[rangeO, H2]
    vOO1 = doo[rangeO, O1]
    vOO2 = doo[rangeO, O2]
    hbond_mat[rangeO, O1] += \
     (np.sum(vOH1 * vOO1, axis=1) / np.sqrt(d2OO1*doh2[rangeO,H1]) > cosacut)
    hbond_mat[rangeO, O2] += \
     (np.sum(vOH2 * vOO2, axis=1) / np.sqrt(d2OO2*doh2[rangeO,H2]) > cosacut)
    # assert (hbond_mat <= 2).all()
    hbond_mat[hbond_mat < 2] = 0
    hbond_mat[hbond_mat >= 2] = 1

    hbond_mat_trj.append(hbond_mat)

    if tenth_nframes > 0:
        if iframe % tenth_nframes == 0 or iframe == nframes-1:
            print(f"hbond_mat_trj finished {iframe/nframes*100}%")

hbond_mat_trj = np.array(hbond_mat_trj)
print("calculating correlation function")
hbond_corr = np.zeros(nframes)
for i in range(nwat):
    for j in range(nwat):
        hbond_corr += autocorrelate(hbond_mat_trj[:,i,j])
np.savetxt(argv[3], hbond_corr / hbond_corr[0])
