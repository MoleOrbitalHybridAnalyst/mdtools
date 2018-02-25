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

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("nbo_out", help='gaussian output with nbo turned on')
    args = parser.parse_args()
 
    debye = 0.20819434

    dipole = {}
    is_dipole = False
    with open(args.nbo_out, "r") as fp:
        for line in fp:
            if is_dipole:
                if re.match(".+BD", line):
                    header = "BD-" + make_atom_str(line[14:17])
                    header += "-" + make_atom_str(line[19:22])
                    dipole[header] = string2array(line[55:73])
                elif re.match(".+LP", line):
                    header = "LP-" + make_atom_str(line[14:17])
                    if header in dipole:
                        dipole[header] += string2array(line[55:73])
                    else:
                        dipole[header] = string2array(line[55:73])
                elif re.match(".+CR", line):
                    header = "CR-" + make_atom_str(line[14:17])
                    dipole[header] = string2array(line[55:73])
                if re.match("\s+deloc", line):
                    if header + "-deloc" in dipole:
                        dipole[header + "-deloc"] += string2array(line[55:73])
                    else:
                        dipole[header + "-deloc"] = string2array(line[55:73])
            if re.match("\s+Residual nuclear charge", line):
                dipole["nuclear"] = string2array(line[55:73])
            if re.match("\s+Net dipole moment", line):
                dipole["total-net"] = string2array(line[55:73])
            if re.match("\s+Delocalization correction", line):
                dipole["total-delocalization"] = string2array(line[55:73])
            if re.match("\s+Total dipole moment", line):
                dipole["total"] = string2array(line[55:73])
            if re.match("\s+Dipole moment analysis", line):
                is_dipole = True
            if re.match("\s+Total dipole moment", line):
                is_dipole = False
    
    deloc = get_match(dipole, ".+deloc")
