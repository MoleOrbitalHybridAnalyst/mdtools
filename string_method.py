import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt

# info about pmf
free = np.loadtxt("free.dat", usecols = [0,1,2])
path = np.loadtxt("initial_path.dat", usecols = [0, 1])
xdim = 200
zmax = 8
xscale = 1.0; yscale = 1.0
# computation settings
update_ends = False      # if update the end points of string
nbeads = 200             # number of beads in string
lr = 5e-3                # learning rate
eps = 1e-6               # stop criteria 
maxiter = 50000          # maximum number of iterations
decay = 5000             # decay lr by 2 every this steps


ydim = int(len(free) / xdim)
if abs(ydim  - len(free) / xdim) > 1e-50:
    print("xdim does not match input data")
    exit(1)

mask = free[:,2] > zmax
free[mask, 2] = zmax
free[:,0] *= xscale; free[:,1] *= yscale
xmin = min(free[:, 0])
xmax = max(free[:, 0])
ymin = min(free[:, 1])
ymax = max(free[:, 1])
free = free.reshape(xdim, -1, 3)

dx = free[1,0,0] - free[0,0,0]
dy = free[0,1,1] - free[0,0,1]
dfx, dfy = np.gradient(free[:,:,2], dx, dy)
dfx_ = RectBivariateSpline(free[:,0,0], free[0,:,1], dfx)
dfy_ = RectBivariateSpline(free[:,0,0], free[0,:,1], dfy)
f_ = RectBivariateSpline(free[:,0,0], free[0,:,1], free[:,:,2])

#def grad_nearest(x, y):
#    xindx = int((x - xmin) / dx)
#    yindx = int((y - ymin) / dy)
#    return [dfx[xindx, yindx], dfy[xindx, yindx]]

def grad(x, y):
    return np.array([dfx_(x,y)[0,0], dfy_(x,y)[0,0]])

def reparametrize_path(path, nbeads):
    # reference https://stackoverflow.com/questions/18244305/how-to-redistribute-points-evenly-over-a-curve

    #equally spaced in arclength
    N = np.transpose(np.linspace(0,1,nbeads))
    
    #number of points on the curve
    n = len(path)
    p1 = path[0,:]
    pend = path[-1,:]
    last_segment = np.linalg.norm(np.subtract(p1,pend))
#    epsilon= 10*np.finfo(float).eps
    
#    #IF the two end points are not close enough lets close the curve
#    if last_segment > epsilon*np.linalg.norm(np.amax(abs(path),axis=0)):
#        path = np.vstack((path,p1))
#        nbeads = nbeads + 1
#    else:
#        print('Contour already closed')
    
    pt = np.zeros((nbeads,2))
    
    #Compute the chordal arclength of each segment.
    chordlen = (np.sum(np.diff(path,axis=0)**2,axis=1))**(1/2)
    #Normalize the arclengths to a unit total
    chordlen = chordlen/np.sum(chordlen)
    #cumulative arclength
    cumarc = np.append(0,np.cumsum(chordlen))
    
    tbins= np.digitize(N,cumarc) # bin index in which each N is in
    
    #catch any problems at the ends
    #tbins[np.where(tbins<=0 | (N <= 0))] = 1
    #tbins[np.where(tbins >= n | (N >= 1))] = n - 1
    tbins[(tbins<=0) | (N <= 0)] = 1
    tbins[(tbins >= n) | (N >= 1)] = n - 1
    
    s = np.divide((N - cumarc[tbins]),chordlen[tbins-1])
    pt = path[tbins,:] + np.multiply((path[tbins,:] - path[tbins-1,:]),(np.vstack([s]*2)).T)

    if not update_ends:
        pt[0] = path[0]
        pt[-1] = path[-1]
    
    return pt

def update_path(path, lr):
    dx = np.gradient(path[:,0])
    dy = np.gradient(path[:,1])
    t = np.array([dx, dy]).T

    path_ = path.copy()
    for i in range(len(path) - 1):
        g = grad(path_[i,0], path_[i,1])
        g_n = g - (g @ t[i]) * t[i] / sum(t[i]**2)
        path_[i] -= lr * g_n
    if not update_ends:
        path_[0] = path[0]
        path_[-1] = path[-1]

    return path_

def save(it):
    np.savetxt("path_" + str(it) + ".dat", path)
    energy = [f_(x,y)[0,0] for x, y in path]
    np.savetxt("energy_" + str(it) + ".dat", energy)

it = 0
while it < maxiter:
    path = reparametrize_path(path, nbeads)
    path_ = update_path(path, lr)
    diff2 = np.sum(( path_ - path)**2 )
    path = path_
    it += 1
    if it % decay == 0: lr /= 2
    if it % (int(maxiter / 20)) == 0:
        print(it, diff2)
        save(it)
    if diff2 < eps:
        print(it, diff2)
        save(it)
        break
