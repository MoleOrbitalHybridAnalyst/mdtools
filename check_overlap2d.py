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
    parser.add_argument('-w','--width',
            help='width of bins in each dimension',default="0.025,0.05")
    parser.add_argument('-r','--ratio',
            help='extend the histo range by this ratio of true range',
            default='0.2,0.2')
    parser.add_argument('-d','--drawline',
            help='draw line at this value',default=0.05)
    args = parser.parse_args()

    c1 = int(args.column.split(',')[0])
    c2 = int(args.column.split(',')[1])
    # n1 = int(args.number.split(',')[0])
    # n2 = int(args.number.split(',')[1])
    w1 = float(args.width.split(',')[0])
    w2 = float(args.width.split(',')[1])
    r1 = float(args.ratio.split(',')[0])
    r2 = float(args.ratio.split(',')[1])
    drawline = float(args.drawline)

    plt.figure()

    cmap = plt.cm.get_cmap('hsv', 20)
    with open(args.list_dir,"r") as fp:
        idir = 0
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
            if len(cv1) == 0:
                continue
            cv1 = np.array(cv1)
            cv2 = np.array(cv2)
            min1 = min(cv1); max1 = max(cv1)
            min2 = min(cv2); max2 = max(cv2)
            extend = (max1 - min1) * r1 /2.0
            min1 -= extend; max1 += extend
            extend = (max2 - min2) * r2 /2.0
            min2 -= extend; max2 += extend
            # w1 = (max1 - min1) / n1
            # w2 = (max2 - min2) / n2
            n1 = int((max1 - min1) / w1)
            n2 = int((max2 - min2) / w2)
            print("min1 = %f max1 = %f min2 = %f max2 = %f"%(min1,max1,min2,max2))
            print("n1 = %d n2 = %d"%(n1,n2))
            histo = np.zeros((n1+1)*(n2+1)).reshape(n2+1,n1+1)
            for (x,y) in zip(cv1,cv2):
                histo[int((y-min2)/w2)][int((x-min1)/w1)] += 1
            histo /= len(cv1)
            x = np.arange(min1,max1+w1,w1)[:n1+1] + w1/2.0
            y = np.arange(min2,max2+w2,w2)[:n2+1] + w2/2.0
            cp = plt.contour(x,y,histo,[drawline],colors=[cmap(idir)])
            cp.levels = [line.rstrip()]
            plt.clabel(cp, cp.levels)
            idir = (idir + 1) % 20

    plt.show()
