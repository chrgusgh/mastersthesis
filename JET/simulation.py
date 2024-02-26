import h5py
import numpy as np
import scipy.constants as const
import sys
sys.path.append('/home/christiang/Skrivbord/DREAM/tools/eqget')

import utils as utils
from EQDSK import EQDSK

class TokamakSimulation:
    def __init__(self, SHOT=1, dBB_cold=1e-4, assimilation=0.2, t_TQ=1.75e-4):
        """
        Initializes the TokamakSimulation class with specific simulation parameters.

        Args:
        - SHOT (int): Index for the tokamak shot to simulate.
        - dBB_cold (float): Initial cold magnetic perturbation parameter.
        - assimilation (float): Percentage of the injected material that enters the tokamak chamber.
        - t_TQ (float): Duration of the thermal quench.
        """
        self.SHOT = SHOT
        self.dBB_cold = dBB_cold
        self.assimilation = assimilation
        self.t_TQ = t_TQ

        # Load shot data and equilibrium configuration
        self.fp = '../JETdata/DDB_Runaways_ArBt_scan.h5'
        self.pulses = ['85021', '85445', '85450', '85451', '85453', '85943']
        self.load_shot_data()

        # Initialize parameters based on loaded data
        self.initialize_parameters()
        
        # Update dependent parameters with provided arguments
        self.update_dependent_parameters()

        # Log input parameters
        self.log_simulation_settings()

    def load_shot_data(self):
        """Loads data for the specified tokamak shot and magnetic equilibrium configuration."""
        self.f = h5py.File(self.fp, 'r')
        self.strr = self.pulses[self.SHOT]
        self.p0 = self.f[self.strr]
        self.mag_eq_fn = f'mag_eq_{self.strr}_data.h5'
        self.eqdsk = EQDSK(f'../JETdata/g_JET_ehtr_{self.strr}_t62.3990_62.4010', override_psilim=2e-3)
        self.equil = h5py.File(f'../JETdata/{self.mag_eq_fn}', 'r')['equil']

    def initialize_parameters(self):
        """Initializes simulation parameters based on loaded shot data."""
        self.ptx = self.equil['ptx'][:]
        self.psi_apRp = np.array([arr.item() for arr in self.equil['psi_apRp'][:]])
        self.a = self.ptx[-1,-1]
        self.R0 = self.eqdsk.R0
        self.psi_physical = self.psi_apRp
        self.psi_n = self.psi_apRp / self.psi_apRp[-1]
        self.plasma_volume = 76.55510067624934  # m^3
        self.Ti = 300  # Assumption for impurities, in Kelvin
        self.process_profiles()
        self.Ip0 = np.abs(self.p0['Ip_A'][:])
        self.q = 1
        self.m_e_eV = const.m_e / const.e
        self.T_i = np.median(self.T0[0])
        self.T_f = 100  # Final temperature in eV
        # TODO: streamline this to accomodate for fluid 1e-9 and isotropic 1e-11.
        self.dt0 = 1e-9
        self.dtmax = 1e-5

    def process_profiles(self):
        """Processes radial profiles for electron density and temperature."""
        self.n0r = self.p0['R_Ne_profile_LIDAR_m'][:] - self.R0
        self.n0r_pos = np.where(self.n0r >= 0)
        self.T0r = self.p0['R_Te_profile_ECE_m'][:] - self.R0
        self.T0r_pos = np.where(self.T0r >= 0)
        # Use filtered data
        self.n0r = self.n0r[self.n0r_pos]
        self.T0r = self.T0r[self.T0r_pos]
        self.n0 = self.p0['Ne_profile_LIDAR_mm3'][:][self.n0r_pos]
        self.T0 = self.p0['Te_profile_ECE_eV'][:][self.T0r_pos]

    def update_dependent_parameters(self):
        """Updates simulation parameters that depend on the initial settings."""
        self.D2 = self.p0['D2_Pam3'][:] * self.assimilation
        self.Ar = self.p0['Ar_Pam3'][:] * self.assimilation
        self.n_D2 = utils.impurity_density(self.D2, self.Ti, self.plasma_volume)
        self.n_Ar = utils.impurity_density(self.Ar, self.Ti, self.plasma_volume)
        self.j_par_at_Bmin = np.abs(self.eqdsk.get_Jpar_at_Bmin(self.psi_n))
        self.r = self.eqdsk.get_r(self.psi_n)
        self.tau_TQ = utils.calculate_tau_TQ(self.t_TQ, self.T_i, self.T_f)
        self.v_th = utils.calculate_v_th(self.T_f, self.m_e_eV)
        self.dBB_re = self.dBB_cold  # Optionally set this based on additional logic
        self.Drr = utils.calculate_Drr(self.R0, self.q, self.dBB_re)

    def log_simulation_settings(self):
        """
        Logs the current simulation settings to a specified log file.
        
        Args:
        - log_file_path (str): Path to the log file where simulation settings will be saved.
        """
        # Prepare the variables and their values in a dictionary
        variables_to_log = {
            #"ptx": self.ptx.tolist(),
            "psi_apRp (scalars)": self.psi_apRp.tolist(),
            "a": self.a,
            "R0": self.R0,
            "psi_physical": self.psi_physical.tolist(),
            "psi_n": self.psi_n.tolist(),
            "plasma_volume": self.plasma_volume,
            "Ti": self.Ti,
            "n0r": self.n0r.tolist(),
            "T0r": self.T0r.tolist(),
            "n0": self.n0.tolist(),
            "T0": self.T0.tolist(),
            "Ip0": self.Ip0.tolist(),
            "D2": self.D2.tolist(),
            "Ar": self.Ar.tolist(),
            "n_D2": self.n_D2,
            "n_Ar": self.n_Ar,
            "j_par_at_Bmin": self.j_par_at_Bmin.tolist(),
            "r": self.r.tolist(),
            "q": self.q,
            "m_e_eV": self.m_e_eV,
            "T_i": self.T_i,
            "T_f": self.T_f,
            "tau_TQ": self.tau_TQ,
            "v_th": self.v_th,
            "dBB_cold": self.dBB_cold,
            "dBB_re": self.dBB_re,
            "Drr": self.Drr,
            "SHOT": self.SHOT,
            "assimilation": self.assimilation,
            "t_TQ": self.t_TQ,
        }

        log_file_path = "simulation_settings.log"
        # Write variables and their values to the specified log file
        with open(log_file_path, "w") as log_file:
            for var_name, value in variables_to_log.items():
                log_file.write(f"{var_name} = {value}\n")
