import matplotlib.pyplot as plt
print("pyplot loaded")
import numpy as np
print("numpy loaded")
import re
from sys import argv
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import zoom


#try:
#    with open("README") as fp:
#        for line in fp:
#            if re.match("wham-2d", line):
#                splits = line.split()
#                xdim = int(splits[4])
#except:
#    print("cannot open README")
#print("first dimension:", xdim)

xdim = int(float(argv[1]))

data = np.loadtxt("errors.dat")
ydim = int(len(data) / xdim)

zmax = 2
pmf_data = np.loadtxt("pmf_averaged.dat")
mask = pmf_data[:,2] > 15
data[:,2][mask] = 2 * zmax

plt.rcParams["figure.figsize"] = (24,18)
plt.rcParams["font.size"] = 55
plt.rcParams["font.family"] = "arial"
plt.rcParams["axes.labelsize"] = 55
xlabel = "CEC Position"
ylabel = "Water Wire Connectivity"

fig, ax = plt.subplots()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
CS = ax.contourf(
        data[:,0].reshape(xdim,-1), data[:,1].reshape(xdim,-1),data[:,2].reshape(xdim,-1), 
        levels = np.linspace(0, zmax * 199 / 200, 200), 
#        levels = np.arange(0, zmax, 0.25), 
        cmap = plt.cm.get_cmap('plasma'))
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('Statistical Uncertainty of PMF (kcal/mol)')
CL = ax.contour(
        data[:,0].reshape(xdim,-1), data[:,1].reshape(xdim,-1),data[:,2].reshape(xdim,-1), 
        levels = np.arange(0, zmax, 0.25), 
        linewidths = 1.5, colors = 'k')
ax.clabel(CL, fontsize=30, inline=1, fmt="%1.2f")
#cbar.add_lines(CS)
#ax.clabel(CS, colors='w', fontsize=10)

##ax3d = plt.gcf().add_subplot(111, projection = '3d')
#DP = ax.pcolor(data[:,0].reshape(xdim,-1), data[:,1].reshape(xdim,-1), data[:,2].reshape(xdim,-1),
##        cmap = plt.cm.get_cmap('plasma'))
##DP = ax3d.plot_surface(data[:,0],data[:,1],data[:,2],cmap=plt.cm.get_cmap('plasma'),rstride=1,cstride=1,linewidth=1)
#cbar = fig.colorbar(DP)
#cbar.ax.set_ylabel('Potential of Mean Force (kcal/mol)')
#plt.show()
plt.savefig("errors.png")
