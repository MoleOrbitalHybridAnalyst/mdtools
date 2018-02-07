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

def to_colvar(df, fname):
    with open(fname, "w") as fp:
        header = "#! FIELDS"
        for _ in df.columns:
            header += " " + _
        print(header, file=fp)
        for index, row in df.iterrows():
            [fp.write(" %8f"%row[_]) for _ in df.columns]
            print("", file=fp)

setattr(pd.DataFrame, "to_colvar", to_colvar)
