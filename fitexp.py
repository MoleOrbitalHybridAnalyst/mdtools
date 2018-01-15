#!/usr/bin/python
from ctypes import (CDLL,POINTER,c_char,c_int,c_float,c_double)
import numpy as np
fitexp=CDLL('/home/chhli/packages/block_average-master/libfitexp.so')
fitting=fitexp.fitting
fitting.argtypes=[POINTER(c_char),POINTER(c_float)]
