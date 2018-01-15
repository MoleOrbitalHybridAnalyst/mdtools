import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_colvar(fname):
    with open(fname, "r") as fp:
        firstline = fp.readline()
    if re.match("#!.+FIELDS", firstline) == None:
        raise NameError('No header found in ' + fname)
    names_ = firstline.split()[2:]
    return pd.read_csv(fname,delim_whitespace=True,names=names_, comment ='#' )
