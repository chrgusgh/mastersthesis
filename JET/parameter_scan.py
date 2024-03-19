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

def save_simulation_data(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ):
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
    # Format filename prefix with dBB and assimilation values
    filename_prefix = f"dBB_{simulation.dBB_cold}_assim_{simulation.assimilation}"

    # Define the target directory based on discharge and filename_prefix
    target_dir = f"../JETresults/parameter_scans/{simulation.discharge}/{filename_prefix}"

    # Ensure the target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Path to the simulation_settings.log file in the current directory
    log_file_path = 'simulation_settings.log'
    # New path for the simulation_settings.log file in the target directory
    new_log_file_path = os.path.join(target_dir, 'simulation_settings.log')

    # Move the simulation_settings.log file
    shutil.move(log_file_path, new_log_file_path)

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

#def save_simulation_data(simulation):
#    """
#    Runs the save_parameter_scan_data.sh script, which stores the current simulation data in
#    ../JETresults/parameter_scans/{pulse number}
#    """
#    subprocess.run(['./save_parameter_scan_data.sh', simulation.discharge, simulation.dBB_cold, simulation.assim], check=True)


def run_DREAM_simulation(argv, SHOT, dBB, assimilation, t_TQ):
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
    simulation = TokamakSimulation(SHOT=SHOT, dBB_cold=dBB, assimilation=assimilation, t_TQ=t_TQ)

    args, settings = get_settings(argv, simulation)

    do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)

    I_RE_final, tau_CQ, I_RE, I_Ohm, I_hot, I_tot, t = get_data(do_TQ, do_CQ, simulation)

    save_simulation_data(simulation, I_RE, I_Ohm, I_hot, I_tot, t, tau_CQ)

    subprocess.run(['./cleanup.sh'])

    return I_RE_final, tau_CQ

def perform_parameter_scan(argv, SHOT, dBB_values, assimilation_values, t_TQ):
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
    resume_i = 0
    resume_j = 8

    # Initialize results matrices
    RE_current_results = np.zeros((len(dBB_values), len(assimilation_values)))
    tau_CQ_results = np.zeros_like(RE_current_results)

    # Assuming dBB_values and assimilation_values are defined

    # Start the parameter scan
    for i, dBB in enumerate(dBB_values):
    # Skip i indices before the resume point
        if i < resume_i:
            continue

        for j, assimilation in enumerate(assimilation_values):
            # If we're at the resume_i, skip j indices before the resume point
            if i == resume_i and j < resume_j:
                continue

            # Now run your simulation or process
            print(f"Processing i={i}, j={j} with dBB={dBB}, assimilation={assimilation}")
            time.sleep(1)
            # Run the simulation with the current parameters
            I_RE_final, tau_CQ = run_DREAM_simulation(argv, SHOT, dBB, assimilation, t_TQ)
            RE_current_results[i, j] = I_RE_final
            tau_CQ_results[i, j] = tau_CQ

            # Reset resume_j to 0 after first use to start from the beginning for subsequent i's
            resume_j = 0

    return dBB_values, assimilation_values, RE_current_results, tau_CQ_results

def plot_contour(dBB_values, assimilation_values, RE_current_results):
    """
    Plots a 2D contour of the final RE current as a function of dBB and assimilation.
    """

    X, Y = np.meshgrid(assimilation_values, dBB_values)
    Z = RE_current_results
    plt.figure(figsize=(10, 6))
    contour_plot = plt.contourf(X, Y, Z, levels=np.linspace(Z.min(), Z.max(), 20), cmap='plasma')
    colorbar = plt.colorbar(contour_plot)
    colorbar.set_label('$I_{RE}$ (MA)') 
    plt.title('Final RE Current as a Function of $\\delta B / B$ and Assimilation')
    plt.xlabel('Assimilation (%)')
    plt.ylabel('$\\delta B / B$')
    plt.show()

def main(argv):
    
    # 0 - 5.
    SHOT = 5
    # Define the ranges for dBB and assimilation
    dBB_range = (1e-4, 1e-2)  # From 1e-4 to 1e-2
    assimilation_range = (1e-2, 1) # From 1% to 100%
    dBB_values = np.linspace(dBB_range[0], dBB_range[1], 10)  # Define dBB values range
    assimilation_values = np.linspace(assimilation_range[0], assimilation_range[1], 10)  # Define assimilation values range
    #TODO: Streamline t_TQ?
    t_TQ = 1e-4
    
    # Perform the parameter scan
    dBB_values, assimilation_values, RE_current_results, tau_CQ_results = perform_parameter_scan(argv,SHOT, dBB_values, assimilation_values, t_TQ)

    #np.savetxt('I_RE_final.txt', RE_current_results, fmt='%f')
    #np.savetxt('tau_CQ.txt', tau_CQ_results, fmt='%f')

    #plot_contour(dBB_values, assimilation_values, RE_current_results)
    print()
    print('Scan completed successfully!')
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
