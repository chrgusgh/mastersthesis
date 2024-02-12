#!/usr/bin/env python3

import argparse
import copy
from driver import a, bar, generate_baseline, parse_args, simulate, MODE_FLUID, MODE_ISOTROPIC
from driver import SCAN_NONE, SCAN_NEON, SCAN_NRE, SCAN_CURRENT
# from visualize import disruption_summary

# from scans import doscan

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import k, e as kb, e
import sys

sys.path.append('/home/christiang/Skrivbord/DREAM/tools/eqget')
from EQDSK import EQDSK
import h5py

# JET data (wikipedia)
#a             = 1.25
#R0            = 2.96
#B0            = 3.45
plasma_volume = 100 # m3
Ti            = 300 #  assumption for impurities
kb            = kb # Boltzmann constant

# Pick one
SHOT = 1

fp      = 'JETdata/DDB_Runaways_ArBt_scan.h5'
pulses  = ['85021', '85445', '85450', '85451', '85453', '85943']
f       = h5py.File(fp, 'r')
p0      = f[pulses[SHOT]]
eq      = EQDSK(f'JETdata/g_JET_ehtr_{pulses[SHOT]}_t62.3990_62.4010', override_psilim=2e-3)

n0r     = p0['R_Ne_profile_LIDAR_m'][:] - eq.R0
n0r_pos = np.where(n0r >= 0)
n0r     = n0r[n0r_pos]
T0r     = p0['R_Te_profile_ECE_m'][:] - eq.R0
T0r_pos = np.where(T0r >= 0)
T0r     = T0r[T0r_pos]
n0      = p0['Ne_profile_LIDAR_mm3'][:]
n0      = n0[n0r_pos]
T0      = p0['Te_profile_ECE_eV'][:]
T0      = T0[T0r_pos]
Ip0     = -p0['Ip_A'][:] # positive current
D2      = p0['D2_Pam3'][:]
Ar      = p0['Ar_Pam3'][:]

# Convert mag equlib data to h5
mag_eq_fn = f'mag_eq_{pulses[SHOT]}_data.h5'
#eq.save_LUKE(mag_eq_fn)
# change
def impurity_density(impurity):
    N = impurity / (kb * Ti)
    n_i = N / plasma_volume
    return n_i

n_D2 = impurity_density(D2)
n_Ar = impurity_density(Ar)

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
            {'name': 'Ar', 'Z': 18, 'n': n_Ar},
            {'name': 'D2', 'Z': 1, 'n': n_D2}
        ], # Impurities
        'Ip0' : Ip0,
        #'nre0': 0,
        #'nre0_r': 0,
        'n0': n0,
        'n0r': n0r,
        'j0': lambda r : np.ones(r.shape),
        'j0r': None,
        #'T0': 1e3,
        #'T0': lambda r : 1000*np.exp(-((r/a)/0.4653)**2),
        'T0': T0,
        'T0r': T0r,
        #'T0': DDB['Te_profile_ECE_eV'][:][:,0],
        #'T0r': DDB['R_Te_profile_ECE_m'][:][:,0],
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
        't_sim': 1e-2,
        'dt0': 2e-10,
        'dtmax': 1e-5,
        'runIoniz': args.runfrom <= 1,
        'runTQ':    args.runfrom <= 2,
        'verboseIoniz': (1 in args.verbose),
        'verboseTQ': (2 in args.verbose)
    }

    return args, settings


def run_disruption_simulation(args, settings):
    ds = generate_baseline(equilibrium=None if args.cylindrical else mag_eq_fn,
        runInit=(args.runfrom <= 0), verboseInit=(0 in args.verbose), **settings)

    doMain, _ = simulate(ds, **settings)

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


