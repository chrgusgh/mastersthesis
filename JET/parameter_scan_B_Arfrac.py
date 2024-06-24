import numpy as np
import matplotlib.pyplot as plt
import time
from simulation import TokamakSimulation
from TCV52717 import get_settings, run_disruption_simulation
from utils import calculate_t_CQ
import sys
import subprocess
import os
import shutil
import h5py

sys.path.append('/home/christiang/Skrivbord/DREAM/tools/eqget')
from EQDSK import EQDSK

def write_simulation_data_txt(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ):
    """
    Saves the results of a DREAM simulation to a set of text files, organizing them within a uniquely named directory based on the simulation parameters.

    Parameters:
    - simulation (object): An object containing attributes of the simulation, including dBB, assimilation, and discharge values.
    - I_RE (numpy.ndarray): 1D array of runaway electron currents measured throughout the simulation.
    - I_Ohm (numpy.ndarray): 1D array of ohmic currents measured throughout the simulation.
    - I_tot (numpy.ndarray): 1D array of total currents (sum of I_RE and I_Ohm) throughout the simulation.
    - t (numpy.ndarray): 1D array of time points corresponding to the current measurements.
    - tau_CQ (float): The calculated current quench time.

    Example file structure for a simulation with discharge '12345', dBB '0.001', and assimilation '50':
    ../JETresults/parameter_scans/12345/dBB_0.001_assim_50/I_RE.txt
    ../JETresults/parameter_scans/12345/dBB_0.001_assim_50/I_Ohm.txt
    ...and so on for I_tot.txt, t.txt, and tau_CQ.txt.
    """
    target_dir = '.'
    # Save arrays to text files within the specific target directory
    np.savetxt(os.path.join(target_dir, 'I_RE.txt'), I_RE, header='I_RE', fmt='%f')
    np.savetxt(os.path.join(target_dir, 'I_Ohm.txt'), I_Ohm, header='I_Ohm', fmt='%f')
    np.savetxt(os.path.join(target_dir, 'I_hot.txt'), I_hot, header='I_hot', fmt='%f')
    np.savetxt(os.path.join(target_dir, 'I_tot.txt'), I_tot, header='I_tot', fmt='%f')
    np.savetxt(os.path.join(target_dir, 't.txt'), t, header='t', fmt='%f')

    # Save tau_CQ separately since it's a scalar
    with open(os.path.join(target_dir, 'tau_CQ.txt'), 'w') as f:
        f.write(f"tau_CQ: {tau_CQ}\n")

def get_data(do_TQ, do_CQ, simulation):
    """
    Gets data from DREAM output files

    Returns:
    - I_RE_final: Final runaway current of CQ phase
    - tau_CQ: Current quench time estimate in seconds.
    """
    I_RE_TQ = do_TQ.eqsys.j_re.current()[:]
    I_RE_CQ = do_CQ.eqsys.j_re.current()[1:]
    I_RE = np.append(I_RE_TQ, I_RE_CQ)

    I_Ohm_TQ = do_TQ.eqsys.j_ohm.current()[:]
    I_Ohm_CQ = do_CQ.eqsys.j_ohm.current()[1:]
    I_Ohm = np.append(I_Ohm_TQ, I_Ohm_CQ)

    I_hot_TQ = do_TQ.eqsys.j_hot.current()[:]
    I_hot_CQ = do_CQ.eqsys.j_hot.current()[1:]
    I_hot = np.append(I_hot_TQ, I_hot_CQ)

    I_tot = I_RE + I_Ohm + I_hot

    t_TQ = do_TQ.grid.t[:]
    t_CQ = do_CQ.grid.t[1:] + t_TQ[-1]
    t = np.append(t_TQ, t_CQ)

    tau_CQ = calculate_t_CQ(I_Ohm, simulation.Ip0, t)

    return I_RE[-1], tau_CQ, I_RE, I_Ohm, I_hot, I_tot, t

def save_simulation_data(simulation):
    """
    Runs the save_parameter_scan_data.sh script, which stores the current simulation data in
    ../JETresults/parameter_scans/{pulse number}
    """
    subprocess.run(['./save_parameter_scan_data.sh', simulation.discharge, str(simulation.B_factor), str(simulation.Arfrac)], check=True)


def run_DREAM_simulation(argv, SHOT, Arfrac, mag_eq_fn, B_factor):
    """
    Executes a DREAM simulation with specified parameters for magnetic perturbation (dBB) and assimilation rate, then processes and saves the simulation results.

    Parameters:
    - argv: Command-line arguments passed to the DREAM simulation, used for configuring the simulation environment.
    - dBB (float): The magnetic perturbation parameter, representing the normalized change in magnetic field strength.
    - assimilation (float): The assimilation rate, indicating the percentage of injected materials that enter the tokamak plasma.

    Returns:
    - I_RE_final (float): The final runaway electron current from the simulation.
    - tau_CQ (float): The calculated current quench duration.

    """
    simulation = TokamakSimulation(SHOT=SHOT, Arfrac=Arfrac, mag_eq_fn=mag_eq_fn, B_factor=B_factor, dBB_cold=2e-3, assimilation=0.1, t_TQ=1e-4)

    args, settings = get_settings(argv, simulation)

    do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)
 
    I_RE_final, tau_CQ, I_RE, I_Ohm, I_hot, I_tot, t = get_data(do_TQ, do_CQ, simulation)

    write_simulation_data_txt(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ)
    
    save_simulation_data(simulation)
    
    subprocess.run(['./cleanup.sh'])

    return 0

def perform_parameter_scan(argv, SHOT, Arfrac_values, B_values):
    """
    Performs a parameter scan over specified ranges of dBB and assimilation values, with an option to resume from a specific point.

    Args:
    - argv: Command line arguments for further configuration.
    - dBB_values: List of dBB values to be scanned.
    - assimilation_values: List of assimilation values to be scanned.
    - t_TQ: Thermal quench time, affecting the simulation.
    - resume_i: Index of dBB_values from which to resume the scan.
    - resume_j: Index of assimilation_values from which to resume the scan.

    Returns:
    - A tuple containing:
        - Array of dBB values.
        - Array of assimilation values.
        - Matrix of final RE current results.
        - Matrix of current quench times.
    """
    
    # Where to start in the simulation
    resume_i = 9
    resume_j = 7

    # Initialize results matrices
    RE_current_results = np.zeros((len(Arfrac_values), len(B_values)))
    tau_CQ_results = np.zeros_like(RE_current_results)

    pulses = ['85021', '85445', '85450', '85451', '85453', '85943']
    discharge = pulses[SHOT]

    eq_fn = f'../JETdata/g_JET_ehtr_{discharge}_t62.3990_62.4010'
    # Load the original magnetic field data
    eq = EQDSK(eq_fn, override_psilim=2e-3)
   
    # Start the parameter scan
    for i, B_factor in enumerate(B_values):
    # Skip i indices before the resume point
        if i < resume_i:
            continue
        
        equil = eq.get_LUKE()
        # Modify magnetic field data
        #B_factor = (i+1) / 3
        equil['ptBx'] = equil['ptBx'] * B_factor
        equil['ptBy'] = equil['ptBy'] * B_factor
        equil['ptBPHI'] = equil['ptBPHI'] * B_factor

        print(B_factor)

        # Save the modified eq object to a new file
        mag_eq_fn = f'../../../../../mnt/DISK4/christiang/data/mag_eq_{discharge}_{B_factor}.h5'

        with h5py.File(mag_eq_fn, 'w') as f:
            f.create_group('equil')

            for key in equil.keys():
                f[f'equil/{key}'] = equil[key]
        
        for j, Arfrac in enumerate(Arfrac_values):
            # If we're at the resume_i, skip j indices before the resume point
            if i == resume_i and j < resume_j:
                continue
            
            # Now run your simulation or process
            print(f"Processing i={i}, j={j} with B={B_factor}, Arfrac={Arfrac}")
            time.sleep(2)
            # Run the simulation with the current parameters
            run_DREAM_simulation(argv, SHOT, Arfrac, mag_eq_fn, B_factor)

            # Reset resume_j to 0 after first use to start from the beginning for subsequent i's
            resume_j = 0

    return 0

def main(argv):

    SHOT = 5
    Arfrac_values = np.linspace(0.001, 0.9999, 21)
    B_values = np.linspace(0.001, 0.9999, 11)
    
    # Perform the parameter scan
    perform_parameter_scan(argv, SHOT, Arfrac_values, B_values)

    print()
    print('Scan completed successfully!')
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
