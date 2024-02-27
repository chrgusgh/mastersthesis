import numpy as np
import matplotlib.pyplot as plt
import time
from simulation import TokamakSimulation
from TCV52717 import get_settings, run_disruption_simulation
from utils import calculate_t_CQ
import sys
import subprocess

def get_data(do_TQ, do_CQ, simulation):
    """
    Gets data from DREAM output files

    Returns:
    - I_RE_final: Final runaway current of CQ phase
    - tau_CQ: Current quench time estimate in seconds.
    """

    I_RE_final = do_CQ.eqsys.j_re.current()[-1]
    
    I_Ohm_TQ = do_TQ.eqsys.j_ohm.current()[:]
    I_Ohm_CQ = do_CQ.eqsys.j_ohm.current()[1:]
    I_Ohm = np.append(I_Ohm_TQ, I_Ohm_CQ)

    #I_p_TQ = do_TQ.eqsys.I_p[:]
    #I_p_CQ = do_CQ.eqsys.I_p[1:]
    #I_p = np.append(I_p_TQ, I_p_CQ)

    t_TQ = do_TQ.grid.t[:]
    t_CQ = do_CQ.grid.t[1:] + t_TQ[-1]
    t = np.append(t_TQ, t_CQ)
    
    tau_CQ = calculate_t_CQ(I_Ohm, simulation.Ip0, t)

    return I_RE_final, tau_CQ

def save_simulation_data(simulation, dBB, assimilation):
    """
    Runs the save_parameter_scan_data.sh script, which stores the current simulation data in
    ../JETresults/parameter_scans/{pulse number}
    """
    dBB_str = str(dBB)
    assimilation_str = str(assimilation)
    subprocess.run(['./save_parameter_scan_data.sh', simulation.discharge, dBB_str, assimilation_str], check=True)

def run_DREAM_simulation(argv, dBB, assimilation):
    """
    Placeholder function for running a DREAM simulation.
    Replace this with actual calls to your DREAM simulation, configured with
    dBB and assimilation parameters, and returning the final RE current.
    """
    simulation = TokamakSimulation(dBB_cold=dBB, assimilation=assimilation)
    args, settings = get_settings(argv, simulation)
    do_TQ, do_CQ = run_disruption_simulation(args, settings, simulation)
    
    time.sleep(1)
    
    I_RE_final, tau_CQ = get_data(do_TQ, do_CQ, simulation)
    save_simulation_data(simulation, dBB, assimilation)

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

    dBB_values = np.linspace(dBB_range[0], dBB_range[1], 20)  # Define dBB values range
    assimilation_values = np.linspace(assimilation_range[0], assimilation_range[1], 100) * 1e-2  # Define assimilation values range
    RE_current_results = np.zeros((len(dBB_values), len(assimilation_values)))  # Initialize results matrix
    tau_CQ_results = RE_current_results
    
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

    np.savetxt('I_RE_final.txt', RE_current_results, fmt='%f')
    np.savetxt('tau_CQ.txt', tau_CQ_results, fmt='%f')

    plot_contour(dBB_values, assimilation_values, RE_current_results)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
