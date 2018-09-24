"""
    read a meta file / diffusion profile
    compute/retrieve position dependent diffusion constants and mftp
"""

import re
import argparse
import numpy as np
from sys import stderr
from scipy.optimize import curve_fit
from numpy import sin,cos,exp

from read_colvar import *
from autocorrelate import *

def parse():

    #parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser( \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--metafile', \
            help = 'input metafile (probably the one used in WHAM)')
    parser.add_argument('--colvar_format', default = 'colvar', \
            help = 'input colvar file format (data or colvar)')
    parser.add_argument('--colvar_name', default = 'cv', \
            help = 'input colvar name')

    parser.add_argument('--dfs_profile', help = 'input/output diffusion file')
    parser.add_argument('--save_acfs', action = 'store_true', \
            help = 'whether to save acfs')
    parser.add_argument('--acf_prefix', default = 'acf', \
            help = 'prefix of acf files')
    parser.add_argument('--acf_ratio', default = 0.01, \
            help = 'ratio of acf to do integration')
    parser.add_argument('--fit_acf', action = 'store_true', \
            help = 'fit acf before doing integration')
    parser.add_argument('--fit_function', \
            default = dft_fit_function, \
            help = 'functional form to be fitted')
    parser.add_argument('--time_name', \
            default = dft_time_name, \
            help = 'variable name of time in fit_function')
    parser.add_argument('--parameter_names', \
            default = dft_parameter_names, \
            help = 'parameter names to be fitted')
    parser.add_argument('--initial_parameters', \
            default = dft_initial_parameters, \
            help = 'parameter names to be fitted')
    parser.add_argument('--save_fitted_acfs', action = 'store_true', \
            help = 'whether to save fitted acfs')
    parser.add_argument('--fitted_acf_prefix', default = 'fitted_acf', \
            help = 'prefix of fitted acf files')

    parser.add_argument('--compute_resistivity', action = 'store_true', \
            help = 'whether to compute resistivity')
    parser.add_argument('--reactant', help = 'cv value of reactant')
    parser.add_argument('--product', help = 'cv value of product')
    parser.add_argument('--npoints_rp', default = 50,
            help = 'number of points to do integration from r to p')

    parser.add_argument('--compute_rate_constant', action = 'store_true', \
            help = 'whether to compute rate constant')
    parser.add_argument('--reactant_bound', \
            help = 'cv value that defines reactant')
    parser.add_argument('--reactant_bound_ts', \
            help = 'cv value that divides reactant and product')
    parser.add_argument('--npoints_rt', default = 50,
            help = 'number of points to do integration from r to ts')

    parser.add_argument('--temp', help = 'temperature')
    parser.add_argument('--kb', default = 0.0019872041, \
            help = 'Boltzmann constant in pmf unit')
    parser.add_argument('--pmf', default = 'pmf.dat', help = 'input pmf')
    parser.add_argument('--pmf_format', default = 'data', \
            help = 'input pmf format (data or colvar)')
    parser.add_argument('--pmf_name', default = 'pmf', \
            help = 'input pmf name')

    parser.add_argument('-v', '--verbose', action = 'store_true', \
            help = 'show info')
    return parser.parse_args()

if __name__ == "__main__":

    dft_fit_function = 'a1*exp(-t**2/b1)+a2*exp(-t**2/b2)'
    dft_time_name = 't'
    dft_parameter_names = 'a1,b1,a2,b2'
    dft_initial_parameters = '0.5,0.01,0.5,0.01'
    args = parse()

    # read meta file
    filenames = []; window_centers = []
    if args.metafile is not None:
        with open(args.metafile) as fp:
            for line in fp:
                # skip empty and comment lines
                if re.match("\s*#", line): continue
                if re.match("\s*$", line): continue
                splits = line.split()
                filenames.append(splits[0])
                window_centers.append(float(splits[1]))
        already_has_dfs = False
    elif args.dfs_profile is not None:
        already_has_dfs = True
        tmp = np.loadtxt(args.dfs_profile)
        if len(tmp.shape) == 1:
            window_centers = [tmp[0]]
            dfs = [tmp[1]]
        else:
            window_centers = tmp[:, 0]
            dfs = tmp[:, 1]
    else:
        raise Exception("No metafile or dfs_profile inputted")

    if not (args.verbose or args.compute_resistivity or \
            args.compute_rate_constant):
        print("Turn on verbose since no resistivity" \
                " or rate constant to be computed")
        args.verbose = True
    if args.save_fitted_acfs and not args.fit_acf:
        args.fit_acf = True
        print("Turn on fit_acf as save_fitted_acfs is demanded")

    # compute diffusion constants
    if args.verbose: print("Diffusion Constants:")
    if not already_has_dfs:
        dfs = []
        for fn, wc in zip(filenames, window_centers):

            # read in colvar
            if args.colvar_format == 'data': 
                colvar_data = np.loadtxt(fn, usecols = [0, 1])
            elif args.colvar_format == 'colvar': 
                colvar_data = read_colvar(fn)[['time', args.colvar_name]].values
            else:
                raise Exception(\
                        "unkown format for colvar files " + args.colvar_format)

            # compute mean and variance
            average = np.mean(colvar_data[:, 1])
            var = np.std(colvar_data[:, 1], ddof = 1)

            # compute acf
            colvar_data[:, 1] -= average
            acf = autocorrelate(colvar_data[:, 1])

            # save acf
            if args.save_acfs:
                np.savetxt(args.acf_prefix + '_' + str(wc), \
                        np.transpose([colvar_data[:, 0], acf]))

            # do fitting
            n = int(float(args.acf_ratio) * len(colvar_data))
            if args.fit_acf:
                if args.fit_function != dft_fit_function:
                    if args.time_name == dft_time_name and \
                       args.parameter_names == dft_parameter_names and \
                       args.initial_parameters == dft_initial_parameters:
                           print('WARNING: you specified fit_function', \
                                 'but did not touch other fit settings', \
                                 file = stderr)
                # define the function to be fitted
                func_string = 'def func(' + args.time_name + ','
                func_string += args.parameter_names
                func_string += '): return ' + args.fit_function
                exec(func_string)

                p0 = np.vectorize(float)(args.initial_parameters.split(','))
                popt, pcov = curve_fit(func, \
                        colvar_data[:n, 0], acf[:n], p0 = p0)

                acf_ = np.array([func(_, *popt) for _ in colvar_data[:, 0]])
                aa = acf_[1:] + acf_[:-1]
                if args.save_fitted_acfs:
                    header = 'functional form:\n'
                    header += '  ' + args.fit_function + '\n'
                    header += 'parameters to be fitted:\n'
                    header += '  ' + args.parameter_names + '\n'
                    header += 'fitted parameters:\n'
                    header += '  ' + ','.join(np.vectorize(str)(popt)) + '\n'
                    header += 'covariance matrix:\n'
                    for row in pcov:
                        header += '  ' + \
                                ','.join(np.vectorize(str)(row)) + '\n'
                    np.savetxt(args.fitted_acf_prefix + '_' + str(wc), \
                            np.transpose([colvar_data[:, 0], acf_]), \
                            header = header)
            else:
                aa = acf[1:] + acf[:-1]

            # compute diffusion constants from acf and var
            dt = colvar_data[1:, 0] - colvar_data[:-1, 0]
            if args.fit_acf:
                # do integ of fitted acf on full range
                integ = sum(dt * aa) / 2.0
            else:
                integ = sum(dt[:n] * aa[:n]) / 2.0
            dfs.append(var**2 / integ)
            if args.verbose:
                print(wc, dfs[-1])
    else:
        for wc, df in zip(window_centers, dfs):
            print(wc, df)

    # save dfs
    if args.dfs_profile is not None and not already_has_dfs:
        np.savetxt(args.dfs_profile, np.transpose([window_centers, dfs]))

    # read temp and pmf and r & p
    if args.compute_resistivity or args.compute_rate_constant:
        kb = float(args.kb)
        if args.temp is None:
            raise Exception("Temperature should be given by --temp")
        kbT = float(args.temp) * kb
        if args.verbose:
            print("kbT =", kbT)
        if args.pmf is None:
            raise Exception("Input pmf file should be given by --pmf")
        if args.pmf_format == 'data':
            pmf_data = np.loadtxt(fn, usecols = [0, 1])
        elif args.pmf_format == 'colvar': 
            pmf_data = read_colvar(fn)[['time', args.pmf_name]].values
        if args.reactant is None:
            raise Exception(
                "CV value of reactant should be given by --reactant")
        r_cv = float(args.reactant)
        if args.product is None:
            raise Exception(
                "CV value of product should be given by --product")
        p_cv = float(args.product)
        npoints_rp = int(args.npoints_rp)

    # compute resistivity
    if args.compute_resistivity:
        cv_range = np.linspace(r_cv, p_cv, num = npoints_rp)
        dcv = cv_range[1:] - cv_range[:-1]
        pmf_interp = np.interp(cv_range, pmf_data[:, 0], pmf_data[:, 1])
        dfs_interp = np.interp(cv_range, window_centers, dfs)
        integrand = np.exp(- pmf_interp / kbT) / dfs_interp
        ii = integrand[1:] + integrand[:-1]
        integ = (ii * dcv) / 2.0
        print("Resistivity:")
        print(integ)

    # compute rate constant
    if args.compute_rate_constant:
        npoints_rt = int(args.npoints_rt)
