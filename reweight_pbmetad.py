# tired of OOP, just want to write functions
from read_colvar import *
import argparse
from sys import stderr
from os import path
import constants as const
import matplotlib.pyplot as plt

def parse():
    parser = argparse.ArgumentParser( \
          formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('colvar', help = 'input COLVAR')
    parser.add_argument('-ph', '--prefix_hills', default = 'HILLS.', \
            help = 'prefix of HILLS files')
    parser.add_argument('-bcv', '--bias_cv_names', 
            help = 'names of CV being biased')
    parser.add_argument('-rcv', '--reweight_cv_names', 
            help = 'names of CV whose PMF to be computed')
    parser.add_argument('-o', '--output',
            help = 'name of output PMF file')
    parser.add_argument('-m', '--method',
            help = 'static/iterative')
    parser.add_argument('-T', '--temperature',
            help = 'temperature')
    parser.add_argument('-prd', '--periodicity', default="None",
            help = 'perodicity of each biased cv')
    parser.add_argument('-nb', '--nbins', default="200",
            help = 'number of bins used when histogramming')
    parser.add_argument('-reg', '--regularization', default="1e-200",
            help = 'substitute this value of zeros in probabilities')
    parser.add_argument('-eps', '--epsilon', default="1e-10",
            help = 'converge criterion for iterations')

    args = parser.parse_args()
    return args

def read_time_series(df, fields):
    assert 'time' == fields[0]
    min_time = df['time'].values[0]
    max_time = df['time'].values[-1]
    dt = df['time'].values[1] - min_time
    trjlen = int((max_time - min_time)/dt)
    values = np.ones([trjlen, len(fields)]) * np.inf
    min_time = df['time'].values[0]
    for l in df[fields].values:
        indx = int((l[0]-min_time)/dt+0.5)
        if indx >= trjlen:
            break
        values[indx] = l

    where_is_missing = np.where(values == np.inf)[0]
    where_is_present = np.where(values != np.inf)[0]

    missing_times =  where_is_missing * dt + min_time
    for time in missing_times:
        print("WARNING: missing time", time, file = stderr)

    # guess the missing values
    if len(where_is_missing):
        for i in range(len(fields)):
            values[where_is_missing,i] = \
                np.interp( \
                    missing_times, 
                    np.arange(
                        min_time, 
                        max_time + dt, dt)[where_is_present],
                    values[where_is_present,i])

    return values

def add_one_gaussian(grid, pos, sigma, height, prd):
    """
    add one gaussian to the existing pmf
    """
    s2 = 2 * sigma**2
    grid[:,1] -= height * np.exp(-(grid[:,0]-pos)**2/s2)
    if prd:
        grid[:,1] -= height * np.exp(-(grid[:,0]-pos-prd)**2/s2)
        grid[:,1] -= height * np.exp(-(grid[:,0]-pos+prd)**2/s2)

def compute_weight(pmfs):
    """
    weights for reweighting
    w = exp(beta*V_PB) = 1 / sum_i exp(-beta*V_G(s_i))
    V_G(s_i) = -F(s_i)
    """
    maxvalue = max([max(_[:,1]) for _ in pmfs])
    for i, f in enumerate(pmfs):
        if i == 0:
            expVPB = np.exp(beta*(f[:,1]-maxvalue))
        else:
            expVPB = (expVPB[...,None] + np.exp(beta*(f[:,1]-maxvalue))[None,:])
    return np.sum(expVPB) / expVPB

def reweight_static(df, dfhs, bcvs, rcvs):
    bcv_ts = read_time_series(df, ['time']+bcvs)

    # determine the ranges of sampled biased cvs
    bmins = np.min(bcv_ts, axis=0)[1:]
    bmaxs = np.max(bcv_ts, axis=0)[1:]

    # 1D pmf of bcvs to be computed from sum hills
    pmfs1d = list()  
    for ib, (b, dfh) in enumerate(zip(bcvs, dfhs)):
        hills = read_time_series(dfh, ['time',b,'sigma_'+b,'height'])
        pmf1d = [np.linspace(bmins[ib],bmaxs[ib],nbins)]
        pmf1d.append(np.zeros(nbins))
        pmf1d = np.array(pmf1d).T
        for h in hills:
            add_one_gaussian(pmf1d, h[1], h[2], h[3], prds[ib])
        pmfs1d.append(pmf1d)

    # this computes the PB weights on the bcv N-dimensional grids
    nd_weights = compute_weight(pmfs1d)
    # compute the weight of each data point in COLVAR
    digits = [ np.digitize(bcv_ts[:,ib+1], pmfs1d[ib][:,0]) \
               for ib in range(len(bcvs)) ]
    for d in digits:
        d[d == nbins] = nbins - 1
    weights = nd_weights[tuple(digits)]

    rcv_ts = read_time_series(df, ['time']+rcvs)
    H, edges = np.histogramdd(rcv_ts[:,1:], bins=nbins, weights=weights)
    mg = np.meshgrid(*[ (e[1:]+e[:-1])/2 for e in edges ], indexing='ij')
    mg.append(H)
    return np.array(mg).reshape(len(rcvs)+1,-1).T

def reweight_iterative(df, dfhs, bcvs, eps, initprob):
    prob = initprob[:,-1].reshape([nbins]*len(bcvs))

    bcv_ts = read_time_series(df, ['time']+bcvs)

    # determine the ranges of sampled biased cvs
    bmins = np.min(bcv_ts, axis=0)[1:]
    bmaxs = np.max(bcv_ts, axis=0)[1:]

    # NOTE I will assume N-d grid of initprob aligns with bins1d
    bins1d = [np.linspace(l,h,nbins) for l,h in zip(bmins,bmaxs)]
    digits = [ np.digitize(bcv_ts[:,ib+1], bins1d[ib]) \
               for ib in range(len(bcvs)) ]
    for d in digits:
        d[d == nbins] = nbins - 1

    hills = [ read_time_series(dfhs[ib], ['time',b,'sigma_'+b,'height']) \
              for ib,b in enumerate(bcvs) ]

    prob_ = np.zeros([nbins]*len(bcvs))
    firsttime = True
    while np.mean((prob - prob_)**2) > eps:
        if not firsttime:
            prob = np.array(prob_)
        pmfs1d = [np.array([b, np.zeros(nbins)]).T for b in bins1d]
        nd_weights = np.ones([nbins]*len(bcvs))
        prob_ = np.zeros([nbins]*len(bcvs))
        ihill = 0
        time_cap = hills[0][0][0]  # time of the first hill
        for istep in range(len(bcv_ts)):

            if ihill == len(hills[0]): break

            # do the weighted histogram
            digits1d = [d[istep] for d in digits] # nbcvs * 1 int
            weight = nd_weights[tuple(digits1d)]  # single float
            prob_[tuple(digits1d)] += weight

            # if hit time_cap, add guassian and continue for looping bcv_ts
            time = bcv_ts[istep][0]
            if time >= time_cap:
                hilllines = [h[ihill] for h in hills] # nbcvs * one line
                time_cap = hilllines[0][0]
                for ib in range(1, len(bcvs)):
                    # make sure times in HILLS files are consistent
                    assert time_cap == hilllines[ib][0]

                for ib in range(len(bcvs)):
                    h = hilllines[ib]
                    add_one_gaussian(pmfs1d[ib], h[1], h[2], h[3], prds[ib])
                expVPB = 1 / compute_weight(pmfs1d)  # exp(-beta*V_PB)
                # reweighting weights =
                # exp(beta*V_PB) * 
                # \int ds exp(-beta*(F(s)+V_PB(s))) / \int ds exp(-beta*F(s))
                nd_weights = np.sum(prob*expVPB) / np.sum(prob) / expVPB
                ihill += 1
        prob_ = prob_ / np.sum(prob_)
        firsttime = False

    mg = np.meshgrid(*bins1d, indexing='ij')
    mg.append(prob_)
    return np.array(mg).reshape(len(bcvs)+1,-1).T


if __name__ == "__main__":

    global beta
    global prds
    global nbins

    args = parse()
    bcvs = args.bias_cv_names.split(',')
    assert len(bcvs) > 1, "why uinsg PBMetaD?"
    rcvs = args.reweight_cv_names.split(',')
    ph = args.prefix_hills
    method = args.method
    nbins = int(args.nbins)
    temp = float(args.temperature)
    beta = 1.0 / const.kb / temp
    reg = float(args.regularization)
    eps = float(args.epsilon)
    assert method in ['static', 'iterative']
    if args.periodicity == "None":
        prds = [None] * len(bcvs)
    else:
        prds_ = args.periodicity.split(',')
        assert len(prds_) == len(bcvs)
        prds = list()
        for _ in prds_:
            try:
                prds.append(float(_))
            except:
                prds.append(eval(_))
    print(f"biased CVs: {bcvs} with periodicity {prds}")
    print(f"CVs to be reweighted: {rcvs}")

#    print(stderr, "WARNING: use with care when CVs have periodicties")

    df_hills = list()
    for b in bcvs:
        df_hills.append(read_colvar(ph+b))

    df_colvar = read_colvar(args.colvar)

    if method == 'static':
        pmf = reweight_static(df_colvar, df_hills, bcvs, rcvs)
    elif method == 'iterative':
        if bcvs != rcvs:
            # TODO do the bcvs == rcvs case first and do an integration
            # to get the PMF of rcvs
            raise NotImplemented
        pmf = reweight_iterative(df_colvar, df_hills, bcvs, eps,
                reweight_static(df_colvar, df_hills, bcvs, bcvs) )

    pmf[:, -1][pmf[:, -1] == 0] = reg
    pmf[:, -1] = -1/beta*np.log(pmf[:, -1])
    np.savetxt(args.output, pmf)
