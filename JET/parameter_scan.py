import numpy as np
import matplotlib.pyplot as plt
import time
from simulation import TokamakSimulation
from TCV52717 import get_settings, run_disruption_simulation
from utils import calculate_t_CQ
import sys
import subprocess
import os

def save_simulation_data(simulation, I_RE, I_Ohm, I_tot, t, tau_CQ):
    """
    Saves the results of a DREAM simulation to a set of text files, organizing them within a uniquely named directory based on the simulation parameters.

    This function creates a directory structure that first categorizes the data by the discharge value of the simulation, and then further organizes it into subdirectories named according to the specific dBB and assimilation values used in the simulation. Each subdirectory contains text files for the runaway electron currents (I_RE), ohmic currents (I_Ohm), total currents (I_tot), time points (t), and the calculated current quench time (tau_CQ).

    Parameters:
    - simulation (object): An object containing attributes of the simulation, including dBB, assimilation, and discharge values.
    - I_RE (numpy.ndarray): 1D array of runaway electron currents measured throughout the simulation.
    - I_Ohm (numpy.ndarray): 1D array of ohmic currents measured throughout the simulation.
    - I_tot (numpy.ndarray): 1D array of total currents (sum of I_RE and I_Ohm) throughout the simulation.
    - t (numpy.ndarray): 1D array of time points corresponding to the current measurements.
    - tau_CQ (float): The calculated current quench time.

    The function dynamically constructs the path to the target directory based on the simulation's discharge attribute and the dBB and assimilation values. It ensures this target directory exists by creating it if necessary. Each parameter array is saved to a separate text file, with the scalar tau_CQ value saved to its own file.

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
    
    # Save arrays to text files within the specific target directory
    np.savetxt(os.path.join(target_dir, 'I_RE.txt'), I_RE, header='I_RE', fmt='%f')
    np.savetxt(os.path.join(target_dir, 'I_Ohm.txt'), I_Ohm, header='I_Ohm', fmt='%f')
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

    I_tot = I_RE + I_Ohm

    t_TQ = do_TQ.grid.t[:]
    t_CQ = do_CQ.grid.t[1:] + t_TQ[-1]
    t = np.append(t_TQ, t_CQ)
    
    tau_CQ = calculate_t_CQ(I_Ohm, simulation.Ip0, t)
    
    return I_RE[-1], tau_CQ, I_RE, I_Ohm, I_tot, t

#def save_simulation_data(simulation):
#    """
#    Runs the save_parameter_scan_data.sh script, which stores the current simulation data in
#    ../JETresults/parameter_scans/{pulse number}
#    """
#    subprocess.run(['./save_parameter_scan_data.sh', simulation.discharge, simulation.dBB_cold, simulation.assim], check=True)


def run_DREAM_simulation(argv, dBB, assimilation):
    """
    Executes a DREAM simulation with specified parameters for magnetic perturbation (dBB) and assimilation rate, then processes and saves the simulation results.

    This function initializes a simulation instance with given dBB and assimilation parameters, retrieves the appropriate settings based on command-line arguments (argv), and runs the disruption simulation. After the simulation, it extracts relevant data including the final runaway electron (RE) current, current quench time (tau_CQ), RE current over time, Ohmic current over time, total current over time, and the corresponding time array. These data are then saved to a structured file, and a cleanup script is executed to prepare for subsequent simulations.

    Parameters:
    - argv: Command-line arguments passed to the DREAM simulation, used for configuring the simulation environment.
    - dBB (float): The magnetic perturbation parameter, representing the normalized change in magnetic field strength.
    - assimilation (float): The assimilation rate, indicating the percentage of injected materials that enter the tokamak plasma.

    The results of the simulation, including the final RE current and the time of the current quench, are returned for further analysis or logging.

    Returns:
    - I_RE_final (float): The final runaway electron current from the simulation.
    - tau_CQ (float): The calculated current quench duration.
    
    """
    simulation = TokamakSimulation(dBB_cold=dBB, assimilation=assimilation)
    
    args, settings = get_settings(argv, simulation)
    
    do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)
    
    time.sleep(1)
    
    I_RE_final, tau_CQ, I_RE, I_Ohm, I_tot, t = get_data(do_TQ, do_CQ, simulation)
    
    save_simulation_data(simulation, I_RE, I_Ohm, I_tot, t, tau_CQ)
    
    subprocess.run(['./cleanup.sh'])

    return I_RE_final, tau_CQ

def perform_parameter_scan(argv, dBB_range, assimilation_range):
    """
    Performs a parameter scan over the specified ranges of dBB and assimilation values.
    
    Args:
    - dBB_range: Tuple containing the min and max values for dBB.
    - assimilation_range: Tuple containing the min and max percentage of assimilation.
    
    Returns:
    - A meshgrid of dBB and assimilation values, and the matrix of RE current results.
    """

    dBB_values = np.linspace(dBB_range[0], dBB_range[1], 10)  # Define dBB values range
    assimilation_values = np.linspace(assimilation_range[0], assimilation_range[1], 50) * 1e-2  # Define assimilation values range
    RE_current_results = np.zeros((len(dBB_values), len(assimilation_values)))  # Initialize results matrix
    tau_CQ_results = RE_current_results
  
    # To deal with a scan that crashes midway
    crash = True
    if crash:
        start_i = 3
        special_start_j_for_i_2 = 22
        # Parameter scan
        for i, dBB in enumerate(dBB_values[start_i:], start=start_i):
        # Determine the starting index for j based on the current i
            if i == 2:
                current_j_start = special_start_j_for_i_2
            else:
                current_j_start = 0  # Reset to start from the first value for other i's

            for j, assimilation in enumerate(assimilation_values[current_j_start:], start=current_j_start):
                I_RE_final, tau_CQ = run_DREAM_simulation(argv, dBB, assimilation)
                RE_current_results[i, j] = I_RE_final
                tau_CQ_results[i, j] = tau_CQ

    else:
        # Parameter scan
        for i, dBB in enumerate(dBB_values):
            for j, assimilation in enumerate(assimilation_values):
                I_RE_final, tau_CQ = run_DREAM_simulation(argv, dBB, assimilation)
                RE_current_results[i, j] = I_RE_final
                tau_CQ_results[i, j] = tau_CQ

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
    # Define the ranges for dBB and assimilation
    dBB_range = (5e-4, 1e-2)  # From 5e-4 to 1e-2
    assimilation_range = (1, 100) # From 1% to 100%

    # Perform the parameter scan
    dBB_values, assimilation_values, RE_current_results, tau_CQ_results = perform_parameter_scan(argv, dBB_range, assimilation_range)

    #np.savetxt('I_RE_final.txt', RE_current_results, fmt='%f')
    #np.savetxt('tau_CQ.txt', tau_CQ_results, fmt='%f')

    plot_contour(dBB_values, assimilation_values, RE_current_results)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
