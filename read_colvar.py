import re
import sys
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

def read_colvar(fname):
    with open(fname, "r") as fp:
        firstline = fp.readline()
    if re.match("#!.+FIELDS", firstline) == None:
        raise RuntimeError('No header found in ' + fname)
    names_ = firstline.split()[2:]
    return pd.read_csv(fname,delim_whitespace=True,names=names_, comment ='#' )

def print_columns(df, stream=sys.stdout):
    for i,c in enumerate(df.columns): print("%d: %s"%(i,c), file=stream)

def print_header(df, stream=sys.stdout):
    header = "#! FIELDS " + " ".join(df.columns)
    print(header, file=stream)

def print_colvar(df, stream=sys.stdout):
    for index, row in df.iterrows():
        [print(" %8f"%row[_], end='', file=stream) for _ in df.columns]
        print("", file=stream)
    
def to_colvar(df, fname):
    with open(fname, "w") as fp:
        print_header(df, fp)
        print_colvar(df, fp)

setattr(pd.DataFrame, "print_columns", print_columns)
setattr(pd.DataFrame, "print_header", print_header)
setattr(pd.DataFrame, "print_colvar", print_colvar)
setattr(pd.DataFrame, "to_colvar", to_colvar)
