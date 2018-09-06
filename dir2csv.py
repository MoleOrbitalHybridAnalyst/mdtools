import argparse
from numpy import argsort, savetxt
from copy import deepcopy

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', help = 'input list of dirs')
    parser.add_argument('--output_csv', help = 'output file name')
    parser.add_argument('--default_cv2', \
            required = True, \
            help = 'default value of cv2 if not present')
    parser.add_argument('--d0', help = 'delimiter of two CVs', default = '-')
    parser.add_argument('--d1', \
            help = 'delimiter of CVs and unrelated info', \
            default = '_')
    parser.add_argument('--orders', \
            help = 'orders of cvs, alphabetic or numeric', \
            default = 'numeric,numeric')
    return parser.parse_args()

def main():
    args = parse()
    orders = args.orders.split(',')

    cv1_cv2 = dict()

    with open(args.input_dir) as fp_in:
        for line in fp_in:
            line = line.rstrip('\n')
            cvs = line.split(args.d1)[0].split(args.d0)
            if cvs[0] in cv1_cv2:
                if len(cvs) == 1:
                    cvs.append(args.default_cv2)
                cv1_cv2[cvs[0]].append(cvs[1])
            else:
                if len(cvs) == 1:
                    cvs.append(args.default_cv2)
                cv1_cv2[cvs[0]] = [cvs[1]]

    keys = list(cv1_cv2.keys())
    if orders[0] == 'alphabetic':
        keys_ = deepcopy(keys)
    elif orders[0] == 'numeric':
        keys_ = []
        for k in keys:
            if k[0] == 'm':
                keys_.append(float('-' + k[1:]))
            else:
                keys_.append(float(k))
    else:
        raise Exception("unkown order type " + orders[0])

    indexes_keys = argsort(keys_)

    all_values = set()
    for k in cv1_cv2:
        values = cv1_cv2[k]
        for v in values:
            all_values.add(v)
    all_values = list(all_values)

    for v in all_values:
        if orders[1] == 'alphabetic':
            all_values_ = deepcopy(all_values)
        elif orders[1] == 'numeric':
            all_values_ = []
            for v in all_values:
                if v[0] == 'm':
                    all_values_.append(float('-' + v[1:]))
                else:
                    all_values_.append(float(v))

    indexes_values = argsort(all_values_)

    results = [['cv2/cv1']]
    for i1 in indexes_values:
        results[0].append(all_values[i1])
    for i0 in indexes_keys:
        result = [keys[i0]]
        values = cv1_cv2[keys[i0]]
        for i1 in indexes_values:
            if all_values[i1] in values: result.append('1')
            else:                        result.append(' ')
        results.append(result)
        #print(result)

    if args.output_csv is None:
        fn_out = args.input_dir + ".csv"
    else:
        fn_out = args.output_csv
    savetxt(fn_out, results, delimiter = ',', fmt = "%s")

if __name__ == "__main__":
    main()
