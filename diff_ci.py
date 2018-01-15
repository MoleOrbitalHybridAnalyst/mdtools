import numpy as np
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', help = 'input ci file 1')
    parser.add_argument('file2', help = 'input ci file 2')
    args = parser.parse_args()

    fp1 = open(args.file1, 'r')
    fp2 = open(args.file2, 'r')
    for (line1, line2) in zip(fp1, fp2):
        s1 = line1.split()
        s2 = line2.split()
        assert s1[0] == s2[0]
        ci1 = {x : float(y) for (x,y) in zip(s1[1::2], s1[2::2])}
        ci2 = {x : float(y) for (x,y) in zip(s2[1::2], s2[2::2])}
        set1 = set(ci1.keys())
        set2 = set(ci2.keys())
        diff_set = set1 - set2
        for _ in diff_set:
            ci2[_] = 0.0
        diff_set = set2 - set1
        for _ in diff_set:
            ci1[_] = 0.0
        d = {x : ci1[x] - ci2[x] for x in ci1}
        print(s1[0], end=' ')
        diff2 = 0.0
        for x in d:
            diff2 += d[x]**2
        print(diff2)
    fp1.close()
    fp2.close()
