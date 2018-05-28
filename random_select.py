# given dir_file 
# generate roster in each $dir

from read_colvar import *
import argparse
from glob import glob

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', help = 'dir file')
    parser.add_argument('-f', '--field', default = 'cv', 
                        help = 'field sampled')
    parser.add_argument('-t', '--timestep', default = 0.001, 
                        help = 'timestep used')
    parser.add_argument('-p', '--prefix', required=True, 
                        help = 'prefix of restart files')
    parser.add_argument('-r', '--range', required=True,
                        help = 'cv range to be selected')
    parser.add_argument('-b', '--nbins', default = 50, 
                        help = 'number of bins')
    parser.add_argument('-s', '--nsamples', default = 1000, 
                        help = 'number of samples')
    parser.add_argument('-e', '--equil_time', default = 50.0, 
                        help = 'equil time in each dir')
    parser.add_argument('-c', '--cv_chk', help = 'checkpoint file for cv')
    
    return parser.parse_args()

def main():
    args = parse()
    timestep = float(args.timestep)
    equil_time = float(args.equil_time)

    # these will be 1D arrays containing each frame's dir, time and cv value:
    dirs = []
    times = []
    cvs = []

    fp_dir = open(args.dir , 'r')
    for d in fp_dir:

        time_white_list = []
        # get all the stepids available in d:
        restart_list = glob(d.strip('\n') + '/' + args.prefix + '*')
        for restart in restart_list:
            time = int(restart.split('.')[-1]) * timestep
            if time >= equil_time: time_white_list.append(time)

        df = read_colvar(d.strip('\n') + '/COLVAR')
        #mask = np.vectorize(select)(df.time)
        mask = df.time.isin(time_white_list)
        dirs.extend(np.repeat(d.strip('\n'), len(df[mask])))
        times.extend(df[mask]['time'].values)
        cvs.extend(df[mask][args.field].values)

    dirs = np.array(dirs)
    times = np.array(times)
    cvs = np.array(cvs)
    
    histo_min, histo_max = np.vectorize(float)(args.range.split(','))
    nbins = int(args.nbins)
    nsamples = int(args.nsamples)
    nsamples_in_one_bin = int(nsamples / nbins)

    selected_indexes = []

    vertexes = np.linspace(histo_min, histo_max, nbins + 1)
    for i in range(nbins):
        low = vertexes[i]; high = vertexes[i+1]
        mask = (cvs >= low) & (cvs < high)

        x = np.arange(len(mask))[mask]
        if len(x) < nsamples_in_one_bin:
            print(
              'WARNING: too few frames (%d) between %f and %f'\
                      %(len(x), low, high),
              file = sys.stderr)
        np.random.shuffle(x)
        x = x[:nsamples_in_one_bin]

        selected_indexes.extend(x)

    fp_dir.seek(0)
    for d in fp_dir:
        mask = dirs[selected_indexes] == d.strip('\n')

        with open(d.strip('\n') + '/roster', 'w') as fp_roster:
            for t in times[selected_indexes][mask]:
                print(args.prefix + str(int(t / args.timestep)), 
                      file = fp_roster)

        if args.cv_chk:
            df = pd.DataFrame(cvs[selected_indexes][mask])
            df.columns = [args.field]
            df.to_colvar(d.strip('\n') + '/' + args.cv_chk)

    fp_dir.close()

if __name__ == '__main__':
    main()
