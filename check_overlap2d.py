import argparse
import matplotlib.pyplot as plt
import numpy as np
import re

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('list_dir', help='a file listing all the dirs')
    parser.add_argument('filename', help='file name to read')
    parser.add_argument('-c','--column',
            help='column numbers of cv1 and cv2',default="1,2")
    parser.add_argument('-n','--number',
            help='number of bins in each dimension',default="10,10")
    parser.add_argument('-r','--ratio',
            help='extend the histo range by this ratio of true range',
            default='0.2,0.2')
    args = parser.parse_args()

    c1 = int(args.column.split(',')[0])
    c2 = int(args.column.split(',')[1])
    n1 = int(args.number.split(',')[0])
    n2 = int(args.number.split(',')[1])
    r1 = float(args.ratio.split(',')[0])
    r2 = float(args.ratio.split(',')[1])

    plt.figure()

    with open(args.list_dir,"r") as fp:
        for line in fp:
            cv1=[]; cv2=[]
            print("checking "+line, end="")
            with open(line.rstrip()+"/"+args.filename,"r") as f:
                for l in f:
                    if re.match("#",l) != None:
                        continue
                    splits = l.split()
                    cv1.append(float(splits[c1]))
                    cv2.append(float(splits[c2]))
            cv1 = np.array(cv1)
            cv2 = np.array(cv2)
            min1 = min(cv1); max1 = max(cv1)
            min2 = min(cv2); max2 = max(cv2)
            extend = (max1 - min1) * r1 /2.0
            min1 -= extend; max1 += extend
            extend = (max2 - min2) * r2 /2.0
            min2 -= extend; max2 += extend
            w1 = (max1 - min1) / n1
            w2 = (max2 - min2) / n2
            histo = np.zeros(n1*n2).reshape(n1,n2)
            for (x,y) in zip(cv1,cv2):
                histo[int((x-min1)/w1)][int((y-min2)/w2)] += 1
            histo /= len(cv1)
            x = np.arange(min1,max1+w1,w1)[:n1] + w1/2.0
            y = np.arange(min2,max2+w2,w2)[:n2] + w2/2.0
            cp = plt.contour(x,y,histo,[0.1])
            cp.levels = [line.rstrip()]
            plt.clabel(cp, cp.levels)

    plt.show()
