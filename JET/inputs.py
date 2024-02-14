#!/usr/bin/env python3

import argparse
import copy
from driver import a, bar, generate_baseline, parse_args, simulate, MODE_FLUID, MODE_ISOTROPIC
from driver import SCAN_NONE, SCAN_NEON, SCAN_NRE, SCAN_CURRENT
# from visualize import disruption_summary

# from scans import doscan

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import k, e as kb, e
import sys

sys.path.append('/home/christiang/Skrivbord/DREAM/tools/eqget')

from EQDSK import EQDSK

# Pick one
SHOT = 1

fp      = '../JETdata/DDB_Runaways_ArBt_scan.h5'
pulses  = ['85021', '85445', '85450', '85451', '85453', '85943']
f       = h5py.File(fp, 'r')
p0      = f[pulses[SHOT]]

eqdsk   = EQDSK(f'../JETdata/g_JET_ehtr_{pulses[SHOT]}_t62.3990_62.4010', override_psilim=2e-3)

# Convert mag equlib data to h5
mag_eq_fn = f'mag_eq_{pulses[SHOT]}_data.h5'
#eq.save_LUKE(mag_eq_fn)
# JET Parameters and more

equil = h5py.File(f'../JETdata/{mag_eq_fn}', 'r')['equil']
ptx = equil['ptx'][:]
psi_apRp = equil['psi_apRp'][:]
# Get each psi_apRp
psi_apRp      = np.array([arr.item() for arr in psi_apRp])
a             = ptx[-1,-1] # 0.901
R0            = eqdsk.R0
psi_physical  = psi_apRp
psi_n         = psi_apRp / psi_apRp[-1]
#B0            = 3.45
plasma_volume = 100 # m3
Ti            = 300 #  assumption for impurities, K
kb            = kb # Boltzmann constant

n0r     = p0['R_Ne_profile_LIDAR_m'][:] - R0
n0r_pos = np.where(n0r >= 0)
n0r     = n0r[n0r_pos]
T0r     = p0['R_Te_profile_ECE_m'][:] - R0
T0r_pos = np.where(T0r >= 0)
T0r     = T0r[T0r_pos]
n0      = p0['Ne_profile_LIDAR_mm3'][:]
n0      = n0[n0r_pos]
T0      = p0['Te_profile_ECE_eV'][:]
T0      = T0[T0r_pos]
Ip0     = -p0['Ip_A'][:] # positive current
D2      = p0['D2_Pam3'][:]
Ar      = p0['Ar_Pam3'][:]

j_par_at_Bmin = eqdsk.get_Jpar_at_Bmin(psi_n)
r             = eqdsk.get_r(psi_n)
#print(r)
#print(j_par_at_Bmin)

def impurity_density(impurity):
    """
    Computes impurity density from the gas law
    
    Input:
    - impurity: in units Pa m^3
    """
    N = impurity / (kb * Ti)
    n_i = N / plasma_volume
    return n_i

n_D2 = impurity_density(D2)
n_Ar = impurity_density(Ar)
