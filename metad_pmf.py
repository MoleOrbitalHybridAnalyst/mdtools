# python3 metad_pmf.py grid.dat gamma T pmf.dat
#   file head in grid.dat should be kept
#   assumes in grid.dat:
#       s1 s2 s3 ... V(s1,s2,s3,...) ...
#   cv numbers are read from FIELDS
#   assume V is in kcal/mol
from scipy.integrate import simps
from sys import argv
import numpy as np
import re

def simps_(xs, y, dim): # \int dxs y(xs)
    integ = y
    for _ in range(dim):
        integ = simps(integ, xs[_])
    return integ

def printhead():
    print("#! FIELDS ",end="",file=fpmf)
    for cvname in cv_names:
        print(cvname,end=" ",file=fpmf)
    print("F",file=fpmf)

if __name__=="__main__":
    num_cv = -1
    V = []; cv_names=[]; cv_mins=[]; cv_maxs=[]; cv_nbins=[]

    with open(argv[1],"r") as fgrid:
        # fetch "FIELDS" to get num_cv cv_names
        for line in fgrid:
            if re.match(".+FIELDS",line)!=None:
                tmp = line.split()
                num_cv = int((len(tmp) - 3)/2)
                for _ in range(2,2+num_cv):
                    cv_names.append(tmp[_])
                break
        assert num_cv != -1, "No FIELDS found in "+argv[1]
        fgrid.seek(0)

        for line in fgrid:
            for cvname in cv_names:
                # fetch "min_cvname" to get cv_mins
                reg=re.match(".+min_"+cvname,line)
                if reg!=None:
                    cv_mins.append(float(line.split()[-1]))
                    continue
                # fetch "max_cvname" to get cv_maxs
                reg=re.match(".+max_"+cvname,line)
                if reg!=None:
                    cv_maxs.append(float(line.split()[-1]))
                    continue
                # fetch "nbins_cvname" to get cv_nbins
                reg=re.match(".+nbins_"+cvname,line)
                if reg!=None:
                    cv_nbins.append(int(line.split()[-1]))
                    continue
        assert len(cv_mins)==num_cv, "min values are missing for some CVs"
        assert len(cv_maxs)==num_cv, "max values are missing for some CVs"
        assert len(cv_nbins)==num_cv, "nbins values are missing for some CVs"
        fgrid.seek(0)
        for icv,cvname in enumerate(cv_names):
            print(cvname,cv_mins[icv],cv_maxs[icv],cv_nbins[icv])

        # construct the V array
        gamma = float(argv[2])
        assert gamma > 1
        temp = float(argv[3])
        dtemp = (gamma - 1) * temp
        kb = 0.0019872041 # boltzmann constant in kcal/mol
        factor = gamma / kb / dtemp;
        for iline,line in enumerate(fgrid):
            if re.match(".+\S+",line)!=None and re.match("^#",line)==None:
                tmp = line.split()
                assert len(tmp) > num_cv, "Cannot read the bias in line "\
                        +str(iline)
                V.append(factor * float(tmp[num_cv]))
        #print(cv_nbins)
        #print(len(V))
        lengthV = len(V)
        V = np.array(V).reshape(cv_nbins)
        fgrid.seek(0)

        # construct cv arrays
        # cvs should have an inverse order of cv
        cvs = []
        for _ in range(num_cv):
            i = num_cv - 1 - _
            cv_array = np.linspace(cv_mins[i],cv_maxs[i],cv_nbins[i])
            cvs.append(cv_array)
        cvs = np.array(cvs)

        # calculate the integration
        shift = kb * temp * np.log(simps_(cvs,np.exp(V),num_cv))
        print("shift =",shift)

        # calculate the shifted V
        V = - kb * temp * V + shift
        V = V.reshape(lengthV)
        iline = 0
        with open(argv[4],"w") as fpmf:
            printhead()
            for line in fgrid:
                if re.match(".+\S+",line)!=None and re.match("^#",line)==None:
                    tmp = line.split()
                    for _ in range(num_cv):
                        print("%18.12f"%float(tmp[_]),end="",file=fpmf)
                    print("%18.12f"%V[iline],file=fpmf)
                    iline += 1
                elif re.match("^#",line)==None or re.match(".+SET",line): \
                        # blank lines or SET lines
                    print(line,end="",file=fpmf)
