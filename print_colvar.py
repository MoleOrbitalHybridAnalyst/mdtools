#! /home/chhli/miniconda3/bin/python3
from read_colvar import *
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('COLVAR', help='input COLVAR')
    parser.add_argument('--fields', help='fields to be printed')
    parser.add_argument('--header', help='also print header', action='store_true')
    return parser.parse_args()

def main():
    args = parse()

    df = read_colvar(args.COLVAR)

    if args.fields is None: 
        df.print_columns()
    else:
        try: df = df[args.fields.split(',')]
        except KeyError as kerr:
            print(kerr, file=sys.stderr)
            print("available fields are:", file=sys.stderr)
            df.print_columns(sys.stderr); exit()

        if args.header:
            df.print_header()
        df.print_colvar()

if __name__=='__main__':
    try:
        main()
    except BrokenPipeError:
        pass
    sys.stderr.close()
