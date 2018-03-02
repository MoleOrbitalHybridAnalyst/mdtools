from read_colvar import *

def ave(s):
    mask1 = pro_cec[2] >= s
    mask2 = pro_cec[2] < s + binsize
    mask = np.logical_and(mask1, mask2)
    return([np.mean(pro_cec[0][mask]), np.mean(pro_cec[1][mask]), s])

if __name__=="__main__":
    df = read_colvar("COLVAR")
    cec_string = "pos_cec"
    cen_string = "pos_cen"
    rot = np.array([[0.92435199,-0.38154081,0.],[0.234192,0.567373,-0.789457],[0.30121,0.729736,0.6138062]])
    binsize = 0.25
    x = df[cec_string + ".x"] - df[cen_string + ".x"]
    y = df[cec_string + ".y"] - df[cen_string + ".y"]
    z = df[cec_string + ".z"] - df[cen_string + ".z"]
    pro_cec = np.dot(rot, [x, y, z])
    zcen = np.arange(min(pro_cec[2]), max(pro_cec[2]), binsize)
    ave_pro = [ave(s) for s in zcen]
    ave_pos = np.dot(ave_pro, rot)
