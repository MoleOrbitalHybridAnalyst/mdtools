from sys import argv, stderr
from os import path
import re

rundir='..'
colvar='COLVAR'
plumed='plumed.dat'
cv1='diste'
cv2=''

if __name__ == '__main__':

   with open(argv[1]) as fp:
      for d in fp:
         cv1at = dict()
         cv1kappa = dict()
         cv2at = dict()
         cv2kappa = dict()
         d = d.rstrip()
         plumed_ = rundir + '/' + d + '/' + plumed
         colvar_ = rundir + '/' + d + '/' + colvar
         if not path.exists(colvar_):
            print(colvar_ + " does not exist", file = stderr)
            continue
         try:
            p = open(plumed_)
         except:
            print("cannot open " + plumed_, file = stderr)

         for line in p:
            if re.match('us:', line):
               args = re.match('.+\sARG=(.+?)\s', line).group(1).split(',')
               ats = re.match('.+\sAT=(.+?)\s', line).group(1).split(',')
               kappas = re.match('.+\sKAPPA=(.+?)\s', line).group(1).split(',')

               cv_at = dict()
               cv_kappa = dict()
               for arg, at, kappa in zip(args, ats, kappas):
                  cv_at[arg] = at
                  cv_kappa[arg] = kappa

               if cv1 in cv_at:
                  if cv1 in cv1at:
                     raise RuntimeError(cv1 + ' appreas on more times in ' + plumed_)
                  cv1at[cv1] = cv_at[cv1]
                  cv1kappa[cv1] = cv_kappa[cv1]
               if cv2 in cv_at:
                  if cv2 in cv2at:
                     raise RuntimeError(cv2 + ' appreas on more times in ' + plumed_)
                  cv2at[cv2] = cv_at[cv2]
                  cv2kappa[cv2] = cv_kappa[cv2]
         p.close()

         if cv1 not in cv1at:
            raise RuntimeError('cannot find ' + cv1 + ' AT in ' + plumed_)
         if cv1 not in cv1kappa:
            raise RuntimeError('cannot find ' + cv1 + ' KAPPA in ' + plumed_)
         if cv2 not in cv2at and cv2 != '':
            raise RuntimeError('cannot find ' + cv2 + ' AT in ' + plumed_)
         if cv2 not in cv2kappa and cv2 != '':
            raise RuntimeError('cannot find ' + cv2 + ' KAPPA in ' + plumed_)
         if cv2 != '':
            print(colvar_, cv1at[cv1], cv2at[cv2], cv1kappa[cv1], cv2kappa[cv2])
         else:
            print(colvar_, cv1at[cv1], cv1kappa[cv1])
