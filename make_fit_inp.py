"""
    read a fit.inp fit.best; output a new fit.inp
"""
import re
import argparse
import numpy as np

def parse():

    parser = argparse.ArgumentParser( \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f',
            help = 'input fit.inp', required = True)
    parser.add_argument('-b',
            help = 'input fit.best', required = True)
    parser.add_argument('-o',
            help = 'output fit.inp', required = True)
    parser.add_argument('-k',
            help = 'field to be fitted', required = True)
    parser.add_argument('-i', default = 'fixed',
            help = 'not fitted if matched')

    return parser.parse_args()

if __name__=="__main__":

    args = parse()

    keys_2b_fitted = args.k.split(",")
    fit_best = np.loadtxt(args.b)

    fp_inp = open(args.f)
    fp_out = open(args.o, 'w')

    indx = 0
    for line in fp_inp:

        if re.match("FIT", line):
            print(line, end = '', file = fp_out)

            # check if this to be fitted
            keywords = line.split()[1:]
            to_be_fitted = False
            if len(keywords) == len(keys_2b_fitted):
                for key, key_ in zip(keywords, keys_2b_fitted):
                    if key != key_:
                        to_be_fitted = False
                        break
                    else:
                        to_be_fitted = True

        elif re.match("\s+$", line):
            print(line, end = '', file = fp_out)
        else:
            if to_be_fitted and not re.match(".+" + args.i, line):
                print(fit_best[indx] / 2, fit_best[indx] * 2, \
                        file = fp_out, end = ' ')
            else:
                print(fit_best[indx], fit_best[indx], \
                        file = fp_out, end = ' ')
            pos_comment = line.index(':')
            print(line[pos_comment:], end = '', file = fp_out)
            indx += 1

    fp_inp.close()
    fp_out.close()
