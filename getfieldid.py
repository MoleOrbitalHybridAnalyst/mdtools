import re
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="input plumed COLVAR file")
    parser.add_argument('field_name', help="field to be looked up")
    args = parser.parse_args()

    fp = open(args.input_file, 'r')
    for line in fp:
        if re.match(".+FIELD", line)!=None:
            splits = line.split()
            for i,_ in enumerate(splits):
                if _==args.field_name:
                    print(i-1)
                    break
            break
    fp.close()

