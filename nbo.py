import numpy as np
import argparse
import re

def get_match(d, pattern):
    return {k: d[k] for k in d if re.match(pattern, k)}

def make_atom_str(string):
    if string[1] == ' ':
        return string[0] + string[2:]
    elif string[1] == '0':
        return string[0] + "1" + string[1:]
    else:
        return string

def string2array(string):
    return np.array(\
            np.vectorize(float)([string[0:6],string[6:12],string[12:18]]))

class nbo:
    def __init__(self, nbo_out):

        self.debye = 0.20819434

        # things will be built
        self.dipole = {}
        self.charges = []

        is_dipole = False
        is_charge = False
        with open(nbo_out, "r") as fp:
            for line in fp:
                # parse dipole
                if is_dipole:
                    if re.match(".+BD", line):
                        header = "BD-" + make_atom_str(line[14:17])
                        header += "-" + make_atom_str(line[19:22])
                        self.dipole[header] = string2array(line[55:73])
                    elif re.match(".+LP", line):
                        header = "LP-" + make_atom_str(line[14:17])
                        if header in self.dipole:
                            self.dipole[header] += string2array(line[55:73])
                        else:
                            self.dipole[header] = string2array(line[55:73])
                    elif re.match(".+CR", line):
                        header = "CR-" + make_atom_str(line[14:17])
                        self.dipole[header] = string2array(line[55:73])
                    if re.match("\s+deloc", line):
                        if header + "-deloc" in self.dipole:
                            self.dipole[header + "-deloc"] += string2array(line[55:73])
                        else:
                            self.dipole[header + "-deloc"] = string2array(line[55:73])
                if re.match("\s+Residual nuclear charge", line):
                    self.dipole["nuclear"] = string2array(line[55:73])
                if re.match("\s+Net dipole moment", line):
                    self.dipole["total-net"] = string2array(line[55:73])
                if re.match("\s+Delocalization correction", line):
                    self.dipole["total-delocalization"] = string2array(line[55:73])
                if re.match("\s+Total dipole moment", line):
                    self.dipole["total"] = string2array(line[55:73])
                if re.match("\s+Dipole moment analysis", line):
                    is_dipole = True
                if re.match("\s+Total dipole moment", line):
                    is_dipole = False

                # parse charge
                if is_charge:
                    if re.match("\s+\* Total \*", line):
                        is_charge = False
                    if re.match("\s+Atom", line): continue
                    splits = line.split()
                    if len(splits) != 7: continue
                    self.charges.append(float(splits[2]))
                if re.match("\s+Summary of Natural Population Analysis", line):
                    is_charge = True
        
        self.deloc = get_match(self.dipole, ".+deloc")
        self.charges = np.array(self.charges)
