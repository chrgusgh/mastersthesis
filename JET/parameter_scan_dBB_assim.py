import numpy as np
import matplotlib.pyplot as plt
import time
from simulation import TokamakSimulation
from TCV52717 import get_settings, run_disruption_simulation
from utils import calculate_t_CQ
import sys
import subprocess
import os

def write_simulation_data_txt(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ):
    """
    Saves the results of a DREAM simulation to a set of text files, organizing them within a uniquely named directory based on the simulation parameters.
    Parameters:
    - simulation (object): An object containing attributes of the simulation, including dBB, assim, and discharge values.
    - I_RE (numpy.ndarray): 1D array of runaway electron currents measured throughout the simulation.
    - I_Ohm (numpy.ndarray): 1D array of ohmic currents measured throughout the simulation.
    - I_tot (numpy.ndarray): 1D array of total currents (sum of I_RE and I_Ohm) throughout the simulation.
    - t (numpy.ndarray): 1D array of time points corresponding to the current measurements.
    - tau_CQ (float): The calculated current quench time.
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

def save_simulation_data(simulation):
    """
    Runs the save_parameter_scan_data.sh script, which stores the current simulation data in
    ../JETresults/parameter_scans/{pulse number}
    """
    subprocess.run(['./save_parameter_scan_data_dBB_assim.sh', simulation.discharge, str(simulation.dBB_cold), str(simulation.assimilation)], check=True)


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

def run_DREAM_simulation(argv, SHOT, dBB, assim, t_TQ):
    """
    Executes a DREAM simulation with specified parameters for mag pert. dBB and assimilation fraction (assim), then processes and saves the simulation results.
    """
    simulation = TokamakSimulation(SHOT=SHOT, dBB_cold=dBB, assimilation=assim, t_TQ=t_TQ)

    args, settings = get_settings(argv, simulation)

    do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)

    I_RE_final, tau_CQ, I_RE, I_Ohm, I_hot, I_tot, t = get_data(do_TQ, do_CQ, simulation)

    write_simulation_data_txt(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ)

    save_simulation_data(simulation)
    
    subprocess.run(['./cleanup.sh'])

    return I_RE_final, tau_CQ

def perform_parameter_scan(argv, SHOT, dBB_values, assim_values, t_TQ):
    """
    Performs a parameter scan over specified ranges of dBB and assim values, with an option to resume from a specific point.

    Args:
    - argv: Command line arguments for further configuration.
    - dBB_values: List of dBB values to be scanned.
    - assim_values: List of assim values to be scanned.
    - t_TQ: Thermal quench time, affecting the simulation.
    - resume_i: Index of dBB_values from which to resume the scan.
    - resume_j: Index of assim_values from which to resume the scan.

    Returns:
    - A tuple containing:
        - Array of dBB values.
        - Array of assim values.
        - Matrix of final RE current results.
        - Matrix of current quench times.
    """

    # Where to start in the simulation
    resume_i = 3
    resume_j = 7


    # Initialize results matrices
    RE_current_results = np.zeros((len(dBB_values), len(assim_values)))
    tau_CQ_results = np.zeros_like(RE_current_results)

    # Assuming dBB_values and assim_values are defined

    # Start the parameter scan
    for i, dBB in enumerate(dBB_values):
        # Skip i indices before the resume point
        if i < resume_i:
            continue

        for j, assim in enumerate(assim_values):
            # If we're at the resume_i, skip j indices before the resume point
            if i == resume_i and j < resume_j:
                continue

            # Now run your simulation or process
            print(f"Processing i={i}, j={j} with dBB={dBB}, assim={assim}")
            time.sleep(2)
            # Run the simulation with the current parameters
            I_RE_final, tau_CQ = run_DREAM_simulation(argv, SHOT, dBB, assim, t_TQ)
            RE_current_results[i, j] = I_RE_final
            tau_CQ_results[i, j] = tau_CQ

            # Reset resume_j to 0 after first use to start from the beginning for subsequent i's
            resume_j = 0

    return dBB_values, assim_values, RE_current_results, tau_CQ_results

def main(argv):
    SHOT = 5
    #dBB_range = (-4, -2)
    assim_range = (0.01, 1)
    #dBB_values = np.logspace(dBB_range[0], dBB_range[1], 11)
    #dBB_values = np.array([1e-4, 2.5e-4, 5e-4, 7.5e-4, 1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2])
    #dBB_values = np.array([1e-4, 5e-4, 1e-3, 5e-3, 1e-2])
    dBB_values = np.array([2.5e-4, 7.5e-4, 2.5e-3, 7.5e-3])
    assim_values = np.linspace(assim_range[0], assim_range[1], 11)
    t_TQ = 1.5e-4

    dBB_values, assim_values, RE_current_results, tau_CQ_results = perform_parameter_scan(argv, SHOT, dBB_values, assim_values, t_TQ)

    print('Scan completed successfully!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
