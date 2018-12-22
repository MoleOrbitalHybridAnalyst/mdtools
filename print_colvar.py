#! /home/chhli/miniconda3/bin/python3
from read_colvar import *
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('COLVAR', help='input COLVAR')
    parser.add_argument('-f', '--fields', help='fields to be printed')
    parser.add_argument('--noheader', 
            help='do not print header', action='store_true')
    parser.add_argument('-r', '--regex', 
            help='match the columns using regex', action='store_true')
    parser.add_argument('-t', '--time',
            help='start_time,time_interval')
    return parser.parse_args()

def main():
    args = parse()

    try:
        df = read_colvar(args.COLVAR)
    except FileNotFoundError as fnferr:
        print(args.COLVAR + ' does not exist', file=sys.stderr); exit()
    except RuntimeError as rterr:
        print(rterr, file=sys.stderr); exit()

    """
        if --fields is absent, print out all the available fields
        if --regex is absent, directly use the input fields to indexing the df
        if --regex is present, use re to find matching fileds, use unique to 
            remove duplicates, use indexes returned by unique to recover the 
            orginal order of the fields
        if no matching fields, throw a KeyError exception
    """

    if args.fields is None: 
        df.print_columns()
    else:
        # get the fields
        raw_f = args.fields.split(',')
        if not args.regex: result_f = raw_f
        else:
            result_f = \
            np.array(
                [col for f in raw_f for col in df.columns if re.match(f,col)])
            indexes = np.unique(result_f, return_index=True)[1]
            result_f = result_f[sorted(indexes)]

        # update df using the fields
        try: 
            if not len(result_f): 
                raise KeyError("no matching fields with " + args.fields)
            times = df.time
            df = df[result_f]
        except KeyError as kerr:
            print(kerr, file=sys.stderr)
            print("available fields are:", file=sys.stderr)
            df.print_columns(sys.stderr); exit()

        # print df
        try:
            if not args.noheader: df.print_header()
            if args.time is None: df.print_colvar()
            else:
                # print cvs at specific times
                min_time, timestep = \
                    np.vectorize(float)(args.time.split(','))
                max_time = max(times)
                timestep_ = times[1] - times[0]
                length = int(round( (max_time - min_time) / timestep )) + 1
                #print(min_time, max_time, timestep, length, len(result_f))
                ts = np.ones(len(result_f) * length) * np.inf
                ts = ts.reshape(length, -1)
                for time, row in zip(times, df.values):
                    indx = int(round( (time - min_time) / timestep ))
                    if abs(indx * timestep + min_time - time) < 0.5 * timestep_:
                        for col, cv in enumerate(row): 
                            ts[indx][col] = cv
                for row in ts:
                    [print(" %8f"%cv, end='') for cv in row]
                    print("")
        except BrokenPipeError:
            pass

        # suppress the 'exception ignored' message
        sys.stderr.close()

if __name__=='__main__':
    main()
