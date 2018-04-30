from read_colvar import *
from autocorrelate import *

boxsize = 19.3887

def pbc(x):
    n = (x + 0.5 * boxsize) // boxsize
    return x - n * boxsize

def do_spec(string):
    x = df[string + '.x'].values
    y = df[string + '.y'].values
    z = df[string + '.z'].values
    vx = x[2:] - x[:-2]
    vy = y[2:] - y[:-2]
    vz = z[2:] - z[:-2]
    vx = np.vectorize(pbc)(vx)
    vy = np.vectorize(pbc)(vy)
    vz = np.vectorize(pbc)(vz)
    dx = df_dipole.dx.values[:len(x)] 
    dy = df_dipole.dy.values[:len(x)]
    dz = df_dipole.dz.values[:len(x)]
    vdx = dx[2:] - dx[:-2]
    vdy = dy[2:] - dy[:-2]
    vdz = dz[2:] - dz[:-2]
    vdx = np.vectorize(pbc)(vdx)
    vdy = np.vectorize(pbc)(vdy)
    vdz = np.vectorize(pbc)(vdz)
    vc = autocorrelate(vx) + autocorrelate(vy) + autocorrelate(vz)
    vc_cross = \
    correlate(vx, vdx-vx) + correlate(vy, vdy-vy) + correlate(vz, vdz-vz) + \
    correlate(vdx-vx, vx) + correlate(vdy-vy, vy) + correlate(vdz-vz, vz)
    vc_rest = autocorrelate(vdx - vx) + \
              autocorrelate(vdy - vy) + \
              autocorrelate(vdz - vz) 
    n = int(len(vc) * 0.01)
    vf = np.real(np.fft.fft(vc[:n]))
    vf_cross = np.real(np.fft.fft(vc_cross[:n]))
    vf_rest = np.real(np.fft.fft(vc_rest[:n]))
    freq = np.fft.fftfreq(n, d = 0.5)
    wn = 33356.40952 * freq
    return [vc, wn, vf, vf_cross, vf_rest]

df = read_colvar("COLVAR")
df_dipole = read_colvar("dipole.dat")

[c_rmdcec, wn_rmdcec, f_rmdcec, f_cross_rmdcec, f_rest_rmdcec] \
        = do_spec("pos_rmdcec")
[c_qc, wn_qc, f_qc, f_cross_qc, f_rest_qc] = do_spec("pos_qc")
[c_pi, wn_pi, f_pi, f_cross_pi, f_rest_pi] = do_spec("pos_pi")
[c_evbcec, wn_evbcec, f_evbcec, f_cross_evbcec, f_rest_evbcec] \
        = do_spec("pos_evbcec")


df_ = pd.DataFrame([[_,__1,__2,__3] for _,__1,__2,__3 in \
        zip(wn_rmdcec,f_rmdcec,f_cross_rmdcec,f_rest_rmdcec)])
df_.columns = ["wavenumber","abs", "cross", "rest", "abs+cross", "total"]
df_[df_.wavenumber>=0].to_colvar("./rmdcec_abs.dat")
df_ = pd.DataFrame([[_,__1,__2,__3,__1+__2,__1+__2+__3] for _,__1,__2,__3 in zip(wn_qc,f_qc,f_cross_qc,f_rest_qc)])
df_.columns = ["wavenumber","abs", "cross", "rest", "abs+cross", "total"]
df_[df_.wavenumber>=0].to_colvar("./qc_abs.dat")
df_ = pd.DataFrame([[_,__1,__2,__3,__1+__2,__1+__2+__3] for _,__1,__2,__3 in zip(wn_pi,f_pi,f_cross_pi,f_rest_pi)])
df_.columns = ["wavenumber","abs", "cross", "rest", "abs+cross", "total"]
df_[df_.wavenumber>=0].to_colvar("./pi_abs.dat")
df_ = pd.DataFrame([[_,__1,__2,__3,__1+__2,__1+__2+__3] for _,__1,__2,__3 in \
        zip(wn_evbcec,f_evbcec,f_cross_evbcec,f_rest_evbcec)])
df_.columns = ["wavenumber","abs", "cross", "rest", "abs+cross", "total"]
df_[df_.wavenumber>=0].to_colvar("./evbcec_abs.dat")

#mask = (wn_rmdcec >= 600) & (wn_rmdcec <= 4000)
mask = wn_rmdcec >= 0
#plt.figure()
#plt.plot(wn_rmdcec[mask], f_rmdcec[mask])
#plt.plot(wn_qc[mask], f_qc[mask])
#plt.plot(wn_pi[mask], f_pi[mask])
#plt.plot(wn_evbcec[mask], f_evbcec[mask])
