#!/usr/bin/env python3

import argparse
import copy
from driver import bar, generate_baseline, parse_args, simulate, MODE_FLUID, MODE_ISOTROPIC
from driver import SCAN_NONE, SCAN_NEON, SCAN_NRE, SCAN_CURRENT
#import inputs as inputs
from simulation import TokamakSimulation
# from visualize import disruption_summary
# from scans import doscan
import matplotlib.pyplot as plt
import numpy as np
import sys

def get_settings(argv, simulation):
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
            {'name': 'Ar', 'Z': 18, 'n': simulation.n_Ar},
            {'name': 'D2', 'Z': 1, 'n': simulation.n_D2}
        ],
        'Ip0' : simulation.Ip0,
        'nre0': 0,
        'nre0_r': 0,
        'n0': simulation.n0,
        'n0r': simulation.n0r,
        #'j0': lambda r : np.ones(r.shape),
        #'j0r': None,
        'j0': simulation.j_par_at_Bmin,
        'j0r': simulation.r,
        'T0': simulation.T0,
        'T0r': simulation.T0r,
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
        't_sim': simulation.t_TQ,
        'dt0': simulation.dt0,
        'dtmax': simulation.dtmax,
        'runIoniz': args.runfrom <= 1,
        'runTQ':    args.runfrom <= 2,
        'verboseIoniz': (1 in args.verbose),
        'verboseTQ': (2 in args.verbose)
    }

    return args, settings


def run_disruption_simulation(args, settings, simulation):
    ds = generate_baseline(equilibrium=None if args.cylindrical else f'../JETdata/{simulation.mag_eq_fn}', simulation=simulation, runInit=(args.runfrom <= 0), verboseInit=(0 in args.verbose), **settings)

    do_TQ, _, do_CQ, _ = simulate(ds1=ds, simulation=simulation, **settings)

    return do_TQ, do_CQ


def main(argv):
    # TODO: streamline this.
    simulation = TokamakSimulation(SHOT=5, dBB_cold=2e-3, assimilation=0.2, t_TQ=1e-4)
    args, settings = get_settings(argv, simulation)

    if args.scan == SCAN_NONE:
        print("Neon density: {} x 10^19 m^-3".format(settings['impurities'][0]['n']/1e19))
        do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)
        #disruption_summary(do_TQ)
    else:
        doscan(args.scan, args, settings, run_disruption_simulation)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
