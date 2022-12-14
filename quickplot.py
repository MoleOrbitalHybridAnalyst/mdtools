def parse():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', 
            nargs='*',
            help='input data files')
    parser.add_argument('-p',  '--plot',
            nargs='*',
            default=None,
            help='what to plot')
    parser.add_argument('-x',  '--xlabel',
            default=None,
            help='x label of plot')
    parser.add_argument('-y',  '--ylabel',
            default=None,
            help='y label of plot')
    parser.add_argument('-s',  '--save',
            default=None,
            help='image file to save')
    parser.add_argument('-l',  '--label',
            nargs='*',
            default=None,
            help='label to each plot')
    parser.add_argument('-g',  '--grid',
            action='store_true',
            help='show grid on plot')
    return parser.parse_args()

def translate(pstr):
    '''
    translate a plot string into data indexes
    e.g. 
        translate('x0:z1') = ("data[0][:,0]", "data[1][:,2]")
        translate('x0:y0-y1') = ("data[0][:,0]", "data[0][:,1]-data[1][:,1]")
    '''
    import re

    p1, p2 = pstr.split(':')

    def sub(s):
        '''
        substitute "x
        '''
        s = re.sub('x([0-9]+)', 'data[\\1][:,0]', s)
        s = re.sub('y([0-9]+)', 'data[\\1][:,1]', s)
        s = re.sub('z([0-9]+)', 'data[\\1][:,2]', s)
        s = re.sub('a([0-9]+)', 'data[\\1][:,3]', s)
        s = re.sub('b([0-9]+)', 'data[\\1][:,4]', s)
        s = re.sub('c([0-9]+)', 'data[\\1][:,5]', s)
        return s
    
    return sub(p1), sub(p2)

if __name__ == "__main__":
    from colors import *
    from numpy import *

    args = parse()

    if args.save is None:
        set_plt_scale(0.4)

    data = [loadtxt(fn) for fn in args.files]

    if args.plot is None:
        args.plot = list()
        for i in range(len(data)):
            args.plot.append(f'x{i}:y{i}')
    if args.label is None:
        args.label = [None] * len(args.plot)
        plt_legend = False
    else:
        plt_legend = True
    
    for pstr, label in zip(args.plot, args.label):
        p1, p2 = translate(pstr)
        if label is None:
            plt.plot(eval(p1), eval(p2))
        else:
            plt.plot(eval(p1), eval(p2), label=label)

    if args.xlabel is not None:
        plt.xlabel(args.xlabel)
    if args.ylabel is not None:
        plt.ylabel(args.ylabel)
    if plt_legend:
        plt.legend(frameon=False)
    if args.grid:
        plt.grid()

    if args.save is None:
        plt.show()
    else:
        plt.savefig(args.save, bbox_inches='tight')
