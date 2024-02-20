#!/usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import sys
sys.path.append('/home/christiang/Skrivbord/DREAM/tools/eqget')

import utils as utils
from EQDSK import EQDSK

# Pick one shot from the list below
SHOT          = 1
fp            = '../JETdata/DDB_Runaways_ArBt_scan.h5'
pulses        = ['85021', '85445', '85450', '85451', '85453', '85943']
f             = h5py.File(fp, 'r')
p0            = f[pulses[SHOT]]
eqdsk         = EQDSK(f'../JETdata/g_JET_ehtr_{pulses[SHOT]}_t62.3990_62.4010', override_psilim=2e-3)
# Convert mag equlib data to h5
mag_eq_fn     = f'mag_eq_{pulses[SHOT]}_data.h5'
#eq.save_LUKE(mag_eq_fn)
equil         = h5py.File(f'../JETdata/{mag_eq_fn}', 'r')['equil']
# JET Parameters and more
ptx           = equil['ptx'][:]
psi_apRp      = equil['psi_apRp'][:]
psi_apRp      = np.array([arr.item() for arr in psi_apRp]) # Get each psi_apRp as a scalar rather than as an arr[scalar]
a             = ptx[-1,-1] # 0.901
R0            = eqdsk.R0
psi_physical  = psi_apRp
psi_n         = psi_apRp / psi_apRp[-1]
#B0            = 3.45
plasma_volume = 76.55510067624934 # m3
Ti            = 300 #  assumption for impurities, K
n0r           = p0['R_Ne_profile_LIDAR_m'][:] - R0
n0r_pos       = np.where(n0r >= 0)
n0r           = n0r[n0r_pos]
T0r           = p0['R_Te_profile_ECE_m'][:] - R0
T0r_pos       = np.where(T0r >= 0)
T0r           = T0r[T0r_pos]
n0            = p0['Ne_profile_LIDAR_mm3'][:]
n0            = n0[n0r_pos]
T0            = p0['Te_profile_ECE_eV'][:]
T0            = T0[T0r_pos]
Ip0           = -p0['Ip_A'][:] # positive current
assimilation  = 0.1 # How much of the injected materials enters the tokamak chamber
D2            = p0['D2_Pam3'][:] * assimilation
Ar            = p0['Ar_Pam3'][:] * assimilation
n_D2          = utils.impurity_density(D2, Ti, plasma_volume)
n_Ar          = utils.impurity_density(Ar, Ti, plasma_volume)
j_par_at_Bmin = eqdsk.get_Jpar_at_Bmin(psi_n)
r             = eqdsk.get_r(psi_n)
# Thermal Quench magnetic perturbation
q             = 1
m_e_eV        = const.m_e / const.e
T_i           = np.median(T0[0]) # initial median plasma temperature in K
T_f           = 100 # eV
t_TQ          = 2e-4 # s, assumption based off of paper by Insulander Björk
tau_TQ        = utils.calculate_tau_TQ(t_TQ, T_i, T_f)
v_th          = utils.calculate_v_th(T_f, m_e_eV)
dBB_cold      = utils.calculate_dBB(a, R0, T_f, t_TQ, v_th)
dBB_re        = utils.calculate_dBB(a, R0, T_f, t_TQ, const.c)
Drr           = utils.calculate_Drr(R0, q, dBB_re)
#dBB           = 1e-3
#Drr           = utils.calculate_Drr(R0, q, dBB)

def log_simulation_settings(log_file_path):
    # Prepare the variables and their values in a dictionary
    variables_to_log = {
        "ptx": ptx,
        "psi_apRp (scalars)": psi_apRp.tolist(),  # Convert np.array to list for JSON-like logging
        "a": a,
        "R0": R0,
        "psi_physical": psi_physical.tolist(),
        "psi_n": psi_n.tolist(),
        "plasma_volume": plasma_volume,
        "Ti": Ti,
        "n0r": n0r.tolist(),
        "T0r": T0r.tolist(),
        "n0": n0.tolist(),
        "T0": T0.tolist(),
        "Ip0": Ip0.tolist(),
        "assimilation": assimilation,
        "D2": D2.tolist(),
        "Ar": Ar.tolist(),
        "n_D2": n_D2,
        "n_Ar": n_Ar,
        "j_par_at_Bmin": j_par_at_Bmin.tolist(),
        "r": r.tolist(),
        "q": q,
        "m_e_eV": m_e_eV,
        "T_i": T_i,
        "T_f": T_f,
        "t_TQ": t_TQ,
        "tau_TQ": tau_TQ,
        "v_th": v_th,
        "dBB_cold": dBB_cold,
        "dBB_re": dBB_re,
        "Drr": Drr
    }

    # Write variables and their values to the specified log file
    with open(log_file_path, "w") as log_file:
        for var_name, value in variables_to_log.items():
            log_file.write(f"{var_name} = {value}\n")

log_file_path = "simulation_settings.log"
log_simulation_settings(log_file_path)
