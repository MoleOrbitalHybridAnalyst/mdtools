import re
import numpy as np
import pandas as pd

def stop_at_special_line(fp, pattern):
    last_pos = fp.tell()
    line = fp.readline()
    isfound = False
    while line:
        if re.match(pattern, line) != None:
            isfound = True; break
        last_pos = fp.tell()
        line = fp.readline()
    fp.seek(last_pos)
    return isfound

class evb_out:
    def __init__(self, evb_out_fn):
        self.fname = evb_out_fn
        self.time_trj = []
        self.states_trj = []
        self.cec_trj = []
        #self.dipole_trj = []
        self.eigenv_trj = []

        fp_evb = open(evb_out_fn, "r")
        while True:
            #parse timestep
            if not stop_at_special_line(fp_evb, "TIMESTEP"): break
            line = fp_evb.readline()
            self.time_trj.append(int(line.split()[1]))
            #parse COMPLEX
            if not stop_at_special_line(fp_evb, "COMPLEX "): break
            self.states_trj.append(evb_states(fp_evb))
            #parse EIGEN_VECTOR
            if not stop_at_special_line(fp_evb, "EIGEN_VECTOR"): break
            line = fp_evb.readline(); line = fp_evb.readline()
            self.eigenv_trj.append(np.vectorize(float)(line.split()))
            #parse CEC
            if not stop_at_special_line(fp_evb, "CEC_COORD"): break
            line = fp_evb.readline(); line = fp_evb.readline()
            self.cec_trj.append(np.vectorize(float)(line.split()))
            ##parse DIPOLE
            #if not stop_at_special_line(fp_evb, "DIPOLE"): break
            #line = fp_evb.readline(); line = fp_evb.readline()
            #self.dipole_trj.append(np.vectorize(float)(line.split()))
            if not line: break
        fp_evb.close()

        self.time_trj = np.array(self.time_trj)
        self.states_trj = np.array(self.states_trj)

class evb_states:
    def __init__(self,fp):
        self.parent = []
        self.shell = []
        self.mol_A = []
        self.mol_B = []
        self.react = []
        self.path = []
        self.extra_cpl = []

        line = fp.readline()
        if line[:7] != "COMPLEX":
            raise Exception("error while parsing COMPLEX")
        #print(line)
        self.nstates = int(re.match(".+:\s+([0-9]+)\s+state", line).group(1))
        line = fp.readline()
        if line[:6] != "STATES":
            raise Exception("error while parsing STATES")
        for _ in range(self.nstates):
            line = fp.readline()
            if re.match("\s+", line) is None:
                break
            splits = line.split()
            if len(splits) != 8:
                raise Exception("incomplete STATES")
            self.parent.append(int(splits[1]))
            self.shell.append(int(splits[2]))
            self.mol_A.append(int(splits[3]))
            self.mol_B.append(int(splits[4]))
            self.react.append(int(splits[5]))
            self.path.append(int(splits[6]))
            self.extra_cpl.append(int(splits[7]))

        self.parent = np.array(self.parent) 
        self.shell = np.array(self.shell)
        self.mol_A = np.array(self.mol_A)
        self.mol_B = np.array(self.mol_B)
        self.react = np.array(self.react)
        self.path = np.array(self.path)
        self.extra_cpl = np.array(self.extra_cpl)

