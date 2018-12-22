# python3 this.py dir [OPTIONS]
# get max block relative entropy (2d) for each window without generating metadata
# for 1d, just save all re.csv and don't need this 
import re
import argparse
import numpy as np
from pathlib import Path
import subprocess

INTERPRETER = 'python3'
BLOCK_RE_SCRIPT = '/home/chhli/share/block_relative_entropy_2d.py'

def parse():

   parser = argparse.ArgumentParser( \
         formatter_class = argparse.ArgumentDefaultsHelpFormatter)

   parser.add_argument('dir', help = 'input list of all windows')

   # inputs for block_relative_entropy
   parser.add_argument('-i', '--input', required = True,
           help = 'input time series')
   parser.add_argument('-f', '--file_format', default = 'colvar', \
           help = 'input colvar file format (data or colvar)')
   parser.add_argument('-n', '--cv_names', default = 'cv,sol', \
           help = 'input colvar names')
   parser.add_argument('-c', '--cv_columns', default = '1,2', \
           help = 'input colvar columns')
   parser.add_argument('-b', '--nblocks', default = 10, \
           help = 'number of blocks the input will be split into')
   parser.add_argument('-w','--widths', \
           help='widths of bins', required = True)
   parser.add_argument('-r','--ratio', \
            help='extend the histo range by this ratio of true range', \
            default = 0.1)

   return parser.parse_args()


if __name__ == '__main__':

   args = parse()

   re_fname = args.dir + "_re.dat"
   path = Path('./' + re_fname)
   if path.is_file():
       print(re_fname + ' exists'); exit();
   fp_re = open(re_fname, 'w')
   print('#! FIELDS window max_re1 max_block1 max_re0 max_block0', file = fp_re)
   window_list = np.genfromtxt(args.dir, dtype = str)
   
   
   for window in window_list:

       print('doing', window, end = ' ')
   
       cv_file = window + '/' + args.input
       path = Path(cv_file)
       if path.is_file(): 
           print('')
       else:
           max_re1 = -1; max_block1 = -1
           max_re0 = -1; max_block0 = -1
           print('not exist, skipped'); continue
       try:
           result = \
            subprocess.check_output( \
                    '%s %s %s -f %s -n %s -c %s -w %s -r %s -b %s'
                    %(INTERPRETER, BLOCK_RE_SCRIPT, cv_file, args.file_format, \
                      args.cv_names, args.cv_columns, args.widths, args.ratio, \
                      args.nblocks) , shell = True)
       except subprocess.CalledProcessError as e:
           print(e); continue;
       result = result.split(b'\n')
       max_re0 = max_re1 = 0.0
       max_block0 = max_block1 = 0
       watershed = int(int(args.nblocks) / 2)
       for line in result:
           if re.match(b'#', line): continue
           if line == b'': continue
           F = line.split()
           if int(F[0]) < watershed:
              if float(F[1]) > max_re0: 
                  max_re0 = float(F[1])
                  max_block0 = int(F[0])
           else:
              if float(F[1]) > max_re1: 
                  max_re1 = float(F[1])
                  max_block1 = int(F[0])

       print(window, max_re1, max_block1, max_re0, max_block0, file = fp_re)

   fp_re.close()
