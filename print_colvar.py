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
    return parser.parse_args()

def main():
    args = parse()

    df = read_colvar(args.COLVAR)

    """
    if --fields is absent, print out all the available fields
    if --regex is absent, directly use the input fields to indexing the df
    if --regex is present, use re to find matching fileds, use unique to remove
        duplicates, use indexes returned by unique to recover the orginal 
        order of the fields
    if no matching fields, throw a KeyError exception
    """

    if args.fields is None: 
        df.print_columns()
    else:
        raw_f = args.fields.split(',')
        if not args.regex: result_f = raw_f
        else:
            result_f = \
            np.array(
                [col for f in raw_f for col in df.columns if re.match(f,col)])
            indexes = np.unique(result_f, return_index=True)[1]
            result_f = result_f[sorted(indexes)]
        try: 
            if not len(result_f): 
                raise KeyError("no matching fields with " + args.fields)
            df = df[result_f]
        except KeyError as kerr:
            print(kerr, file=sys.stderr)
            print("available fields are:", file=sys.stderr)
            df.print_columns(sys.stderr); exit()

        if not args.noheader: df.print_header()
        df.print_colvar()

if __name__=='__main__':
    try:
        main()
    except BrokenPipeError:
        pass
    sys.stderr.close()
