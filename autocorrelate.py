import numpy as np
from scipy import signal
from scipy import fftpack

def autocorrelate(v):
    n = len(v)
    # padding
    if n % 2:
        v_pad = np.zeros(n * 2 - 1)
    else:
        v_pad = np.zeros(n * 2)
    v_pad[n//2 : n//2 + n] = v
    return (signal.fftconvolve(v_pad,v_pad[::-1],mode='same'))[(-n):] / np.arange(n,0,-1)

def correlate(v, s): 
    """
    compute <v(t)s(0)>
    """
    assert len(v) == len(s)
    n = len(v)
    if n % 2:
        v_pad = np.zeros(n * 2  - 1)
        s_pad = np.zeros(n * 2  - 1)
    else:
        v_pad = np.zeros(n * 2)
        s_pad = np.zeros(n * 2)
    v_pad[n//2 : n//2 + n] = v
    s_pad[n//2 : n//2 + n] = s
    return (signal.fftconvolve(v_pad,s_pad[::-1],mode='same'))[(-n):] / np.arange(n,0,-1)

def autocorrelate_vector(v):
    """
    compute \\langle \\vec{v}(0)\\cdot\\vec{v}(t) \\rangle
    """
    result = np.zeros(len(v))
    for d in range(len(v[0])):
        result += autocorrelate(v[:,d])
    return result
