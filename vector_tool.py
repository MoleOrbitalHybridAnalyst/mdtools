import numpy as np
import re

def array2string(arr, sep = ' '):
    return sep.join(np.vectorize(str)(arr))

def string2array(s, dtype = float):
    if re.search(',', s):
        splits = s.split(',')
    else:
        splits = s.split()
    return np.array(np.vectorize(dtype)(splits))
