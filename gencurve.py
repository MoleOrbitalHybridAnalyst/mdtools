from read_colvar import *
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cv_file', help='input colvar file',
            default = 'COLVAR')
    parser.add_argument("--cec_name", help='name of action for the cec',
            default = 'pos_cec')
    parser.add_argument("--cen_name", help='name of action for the cen',
            default = 'pos_cen')
    parser.add_argument("--rotx", help='x vector of the channel',
            default = '0.92435199,-0.38154081,0')
    parser.add_argument("--roty", help='y vector of the channel',
            default = '0.234192,0.567373,-0.789457')
    parser.add_argument("--rotz", help='z vector of the channel',
            default = '0.30121,0.729736,0.6138062')
    parser.add_argument("--binsize", help='binsize for making slices',
            default = '0.25')
    return parser.parse_args()

def ave(s):
    mask1 = pro_cec[2] >= s
    mask2 = pro_cec[2] < s + binsize
    mask = np.logical_and(mask1, mask2)
    return([np.mean(pro_cec[0][mask]), np.mean(pro_cec[1][mask]), s])

if __name__=="__main__":
    args = parse()
    df = read_colvar(args.cv_file)
    cec_string = args.cec_name
    cen_string = args.cen_name
    rot = np.array([
        np.vectorize(float)(args.rotx.split(",")),
        np.vectorize(float)(args.roty.split(",")),
        np.vectorize(float)(args.rotz.split(","))
        ])
    binsize = float(args.binsize)
    x = df[cec_string + ".x"] - df[cen_string + ".x"]
    y = df[cec_string + ".y"] - df[cen_string + ".y"]
    z = df[cec_string + ".z"] - df[cen_string + ".z"]
    pro_cec = np.dot(rot, [x, y, z])
    zcen = np.arange(min(pro_cec[2]), max(pro_cec[2]), binsize)
    ave_pro = [ave(s) for s in zcen]
    ave_pos = np.dot(ave_pro, rot)
