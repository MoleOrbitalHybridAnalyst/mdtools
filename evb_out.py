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

class evb_out(object):
    def __init__(self, evb_out_fn):
        self.fname = evb_out_fn
        self.time_trj = []
        self.ene_environment_trj = []
        self.ene_complex_trj = []
        self.ene_total_trj = []
        self.states_trj = []
        self.diag_trj = []
        self.offdiag_trj = []
        self.extra_cpl_trj = []
        self.cec_trj = []
        #self.dipole_trj = []
        self.eigenv_trj = []

        fp_evb = open(evb_out_fn, "r")
        while True:
            #parse timestep
            if not stop_at_special_line(fp_evb, "TIMESTEP"): break
            line = fp_evb.readline()
            self.time_trj.append(int(line.split()[1]))
            #parse ENE_ENVIRONMENT
            if not stop_at_special_line(fp_evb, "ENE_ENVIRONMENT "): break
            line = fp_evb.readline()
            self.ene_environment_trj.append(float(line.split()[1]))
            #parse ENE_COMPLEX
            if not stop_at_special_line(fp_evb, "ENE_COMPLEX "): break
            line = fp_evb.readline()
            self.ene_complex_trj.append(float(line.split()[1]))
            #parse ENE_TOTAL
            if not stop_at_special_line(fp_evb, "ENE_TOTAL "): break
            line = fp_evb.readline()
            self.ene_total_trj.append(float(line.split()[1]))
            #parse COMPLEX
            if not stop_at_special_line(fp_evb, "COMPLEX "): break
            self.states_trj.append(evb_states(fp_evb))
            #parse DIAGONAL
            if not stop_at_special_line(fp_evb, "DIAGONAL "): break
            self.diag_trj.append(evb_diag(fp_evb,self.states_trj[-1].nstates))
            #parse OFF-DIAGONAL
            if not stop_at_special_line(fp_evb, "OFF-DIAGONAL "): break
            self.offdiag_trj.append(evb_offdiag(fp_evb,self.states_trj[-1].nstates))
            #parse EXTRA-COUPLING
            if not stop_at_special_line(fp_evb, "EXTRA-COUPLING "): break
            self.extra_cpl_trj.append(evb_extra_cpl(fp_evb))
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

    def get_offdiag_graph(self, istep):
        """
            get {(mol_A, mol_B): energy} at istep
        """
        offdiag = self.offdiag_trj[istep]
        extra_cpl = self.extra_cpl_trj[istep]
        states = self.states_trj[istep]
        graph = dict()
        for istate, (mol_A, mol_B) in enumerate(zip(states.mol_A,states.mol_B)):
            if mol_A == -1: continue
            graph[mol_A, mol_B] = offdiag.energy[istate - 1]
        for I, J, energy in zip(extra_cpl.I,extra_cpl.J,extra_cpl.energy):
            graph[states.mol_B[I],states.mol_B[J]] = energy

        return graph

    def get_diag_dict(self, istep):
        """
            get {mol_B: energy} at istep
        """
        diag = self.diag_trj[istep]
        states = self.states_trj[istep]
        d = dict()
        for istate, mol_B in enumerate(states.mol_B):
            d[mol_B] = diag.total[istate]
        return d

    def get_hb_acceptors(self, istep, state_id):
        states = self.states_trj[istep]
        return list(np.arange(len(states.mol_A))[states.parent == state_id])

    def get_hb_donor(self, istep, state_id):
        states = self.states_trj[istep]
        return states.parent[state_id]

class evb_states(object):
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

class evb_diag(object):
    def __init__(self, fp, nstates):
        self.total = []
        self.vdw = []
        self.coul = []
        self.bond = []
        self.angle = []
        self.dihedral = []
        self.improper = []
        self.kspace = []
        self.repulsive = []

        line = fp.readline()
        if line[:8] != "DIAGONAL":
            raise Exception("error while parsing DIAGONAL")
        for _ in range(nstates):
            line = fp.readline()
            if re.match("\s+", line) is None:
                break
            splits = line.split()
            if len(splits) != 10:
                raise Exception("incomplete DIGONAL")
            self.total.append(float(splits[1]))
            self.vdw.append(float(splits[2]))
            self.coul.append(float(splits[3]))
            self.bond.append(float(splits[4]))
            self.angle.append(float(splits[5]))
            self.dihedral.append(float(splits[6]))
            self.improper.append(float(splits[7]))
            self.kspace.append(float(splits[8]))
            self.repulsive.append(float(splits[9]))

        self.total = np.array(self.total)
        self.vdw = np.array(self.vdw)
        self.coul = np.array(self.coul)
        self.bond = np.array(self.bond)
        self.angle = np.array(self.angle)
        self.dihedral = np.array(self.dihedral)
        self.improper = np.array(self.improper)
        self.kspace = np.array(self.kspace)
        self.repulsive = np.array(self.repulsive)

class evb_offdiag(object):
    def __init__(self, fp, nstates):
        self.energy = []
        self.A_Rq = []
        self.Vij_const = []
        self.Vij = []

        line = fp.readline()
        if line[:12] != "OFF-DIAGONAL":
            raise Exception("error while parsing OFF-DIAGONAL")
        for _ in range(1, nstates):
            line = fp.readline()
            if re.match("\s+", line) is None:
                break
            splits = line.split()
            if len(splits) != 5:
                raise Exception("incomplete OFF-DIGONAL")
            self.energy.append(float(splits[1]))
            self.A_Rq.append(float(splits[2]))
            self.Vij_const.append(float(splits[3]))
            self.Vij.append(float(splits[4]))

        self.energy = np.array(self.energy)
        self.A_Rq = np.array(self.A_Rq)
        self.Vij_const = np.array(self.Vij_const)
        self.Vij = np.array(self.Vij)

class evb_extra_cpl(object):
    def __init__(self, fp):
        self.I = []
        self.J = []
        self.energy = []

        line = fp.readline()
        if line[:14] != "EXTRA-COUPLING":
            raise Exception("error while parsing EXTRA-COUPLING")
        self.nextra_cpl = int(re.match("EXTRA-COUPLING\s+([-0-9]+)", line).group(1))
        for _ in range(self.nextra_cpl):
            line = fp.readline()
            if re.match("\s+", line) is None:
                break
            splits = line.split()
            if len(splits) != 3:
                raise Exception("incomplete EXTRA-COUPLING")
            self.I.append(int(splits[0]))
            self.J.append(int(splits[1]))
            self.energy.append(float(splits[2]))

        self.I = np.array(self.I)
        self.J = np.array(self.J)
        self.energy = np.array(self.energy)
