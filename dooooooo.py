from scipy.integrate import simps
from sys import argv
import numpy as np
import re

#def multisimps(xs, y, dim): # \int dxs y(xs)
#    integ = y
#    for _ in range(1,dim):
#        integ = simps(integ, xs[:,_])

def multisimps(xs, y, dim): # \int dxs y(xs)
    integ = y
    for _ in range(dim):
        integ = simps(integ, xs[_])
    return integ

if __name__=="__main__":
    #cvs = []
    #V= []
    ##for cv in np.ndindex((100,100,100)):
    #for cv in np.ndindex((100,100)):
    #    cv = np.array(cv)
    #    cvs.append(cv)
    #    #if (cv[0]-50)**2+(cv[1]-50)**2+(cv[2]-50)**2<=100:
    #    if (cv[0]-30)**2+(cv[1]-30)**2<=100:
    #        V.append(1)
    #    else:
    #        V.append(0)
    #cvs = np.array(cvs)
    #print(cvs[:,0])
    #print(cvs[:,1])
    ##print(multisimps(cvs,V,2))
    #print(multisimps(np.array([range(100),range(100)]),V,2))

    x = np.linspace(0,1,200)
    y = np.linspace(0,1,300)
    z = np.linspace(0,1,400)
    V = []
    for _ in y:
        for __ in x:
            for ___ in z:
                if (_-0.5)**2+(__-0.5)**2+(___-0.5)**2<=0.25:
                    V.append(1)
                else:
                    V.append(0)
    print(len(V))
    #V = np.array(V).reshape((300,200))
    V = np.array(V).reshape((400,300,200))
    #z = np.cos(x)**4 + np.sin(y[:,None])**2
    print(multisimps(np.array([x,y,z]),V,3))

