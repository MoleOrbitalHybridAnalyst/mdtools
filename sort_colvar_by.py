from read_colvar import *
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('colvar', help = 'input COLVAR')
parser.add_argument('-f', help = 'field to be sorted')
parser.add_argument('-o', default = 0, 
        help = 'order descending if zero')

args = parser.parse_args()

df = read_colvar(args.colvar)
df.sort_values(by = args.f, ascending = args.o).to_colvar(args.colvar, dtype = str)
