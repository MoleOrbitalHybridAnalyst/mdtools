from collections import defaultdict
import numpy as np
import argparse
import re

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("nbo_out", help='gaussian output with nbo turned on')
    args = parser.parse_args()

    is_dipole = False
    with open(args.nbo_out, "r") as fp:
        for line in fp:
            # TODO construct the dicts
            if re.match("\s+Diople moment analysis", line) != None:
                is_dipole = True
            if re.match("\s+Total dipole moment", line) != None:
                is_dipole = False
