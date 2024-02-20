#!/usr/bin/env python3

import argparse
import copy
from driver import bar, generate_baseline, parse_args, simulate, MODE_FLUID, MODE_ISOTROPIC
from driver import SCAN_NONE, SCAN_NEON, SCAN_NRE, SCAN_CURRENT
import inputs as inputs
# from visualize import disruption_summary
# from scans import doscan
import matplotlib.pyplot as plt
import numpy as np
import sys

def get_settings(argv):
    args = parse_args(argv)
    ext = args.extension
    PREFIX = 'output/52717'

    if ext:
        ext = '_'+ext

    settings = {
        'mode': args.mode,
        'reltol': 2e-6,
        'prefix': PREFIX,
        'extension': ext,
        'impurities': [
            #{'name': 'Ne', 'Z': 10, 'n': bar(6.2)}
            #{'name': 'Ne', 'Z': 10, 'n': 1e19},
            {'name': 'Ar', 'Z': 18, 'n': inputs.n_Ar},
            {'name': 'D2', 'Z': 1, 'n': inputs.n_D2}
        ],
        'Ip0' : inputs.Ip0,
        'nre0': 0,
        'nre0_r': 0,
        'n0': inputs.n0,
        'n0r': inputs.n0r,
        #'j0': lambda r : np.ones(r.shape),
        #'j0r': None,
        'j0': inputs.j_par_at_Bmin,
        'j0r': inputs.r,
        'T0': inputs.T0,
        'T0r': inputs.T0r,
        'E0': 1,
        'E0r': None,
        'tauwall': None,
        #'tauwall': 0.01,
        #'Vloop': 1,
        'Vloop': None,
        'withfre': args.withfre,
        # Time settings
        #'t_ioniz':  1e-6,
        #'nt_ioniz': 8000,
        #'t_sim': 1e-2,
        't_sim': inputs.tau_TQ,
        'dt0': 2e-13,
        'dtmax': 1e-6,
        'runIoniz': args.runfrom <= 1,
        'runTQ':    args.runfrom <= 2,
        'verboseIoniz': (1 in args.verbose),
        'verboseTQ': (2 in args.verbose)
    }

    return args, settings


def run_disruption_simulation(args, settings):
    ds = generate_baseline(equilibrium=None if args.cylindrical else f'../JETdata/{inputs.mag_eq_fn}',
        runInit=(args.runfrom <= 0), verboseInit=(0 in args.verbose), **settings)

    doMain, _ = simulate(ds1=ds, Drr = inputs.Drr, dBB0 = inputs.dBB_cold, **settings)

    return doMain


def main(argv):
    args, settings = get_settings(argv)

    if args.scan == SCAN_NONE:
        print("Neon density: {} x 10^19 m^-3".format(settings['impurities'][0]['n']/1e19))
        doMain = run_disruption_simulation(args, settings)
        disruption_summary(doMain)
    else:
        doscan(args.scan, args, settings, run_disruption_simulation)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


