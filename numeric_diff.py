import argparse
from read_colvar import *

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', help='input file1')
    parser.add_argument('file2', help='input file2')
    parser.add_argument('--col1', help='columns to be compared in file1',
            default = '0')
    parser.add_argument('--col2', help='columns to be compared in file2',
            default = '0')
    parser.add_argument('--prds', help='periods of the fields')
    parser.add_argument('--result', help='output file of differences')
    return parser.parse_args()

def pbc_diff(x,y,p):
    if args.prds == None:
        return x - y
    else:
        d = x - y
        k = np.floor(d / p + 0.5)
        return d - k * p 

if __name__=='__main__':
    args = parse()
    col1 = np.vectorize(int)(args.col1.split(','))
    col2 = np.vectorize(int)(args.col2.split(','))
    prds = np.zeros(len(col1))
    assert len(col1) == len(col2)
    if args.prds != None:
        prds = np.vectorize(float)(args.prds.split(','))
        assert len(col1) == len(prds)
    if args.file1 == args.file2:
        fp1 = open(args.file1, "r")
    else:
        fp1 = open(args.file1, "r")
        fp2 = open(args.file2, "r")
    result = []
    if args.file1 == args.file2:
        for line1 in fp1:
            if re.match("\s*#", line1): continue
            values1 = np.vectorize(float)(line1.split())[col1]
            values2 = np.vectorize(float)(line1.split())[col2]
            # diff = values1 - values2
            diff = pbc_diff(values1, values2, prds)
            ph = []
            ph.extend(diff)
            ph.append(np.sum(diff**2))
            result.append(ph)
    else:
        for line1,line2 in zip(fp1,fp2):
            if re.match("\s*#", line1): continue
            if re.match("\s*#", line2): continue
            values1 = np.vectorize(float)(line1.split())[col1]
            values2 = np.vectorize(float)(line2.split())[col2]
            # diff = values1 - values2
            diff = pbc_diff(values1, values2, prds)
            ph = []
            ph.extend(diff)
            ph.append(np.sum(diff**2))
            result.append(ph)
    fp1.close()
    if args.file1 != args.file2: fp2.close()
    diff2 = [_[-1] for _ in result]
    print("largest deviation =",np.max(diff2),"at line",np.argmax(diff2)+1)
    print("average deviation =",np.mean(diff2))
    print("stddev of deviation =",np.std(diff2))
    if args.result != None:
        df = pd.DataFrame(np.array(result))
        columns = []
        for col1_,col2_ in zip(col1,col2):
            columns.append("file1."+str(col1_)+"-file2."+str(col2_))
        columns.append("diff2")
        df.columns = columns
        df.to_colvar(args.result)
