from read_colvar import *
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('colvar', help = 'input COLVAR')
parser.add_argument('-f', help = 'field to be filtered')
parser.add_argument('-c', help = 'condition for keeping')

args = parser.parse_args()

df = read_colvar(args.colvar)
exec("mask = (" + args.c%"df[args.f]" + ")")
df[mask].to_colvar(args.colvar)
