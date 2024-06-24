import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from scipy.ndimage import uniform_filter1d

def moving_average(data, window_size):
    """Applies a moving average filter to smooth the data."""
    return uniform_filter1d(data, size=window_size, mode='nearest')

def calculate_derivative(data, time):
    """Calculates the derivative of data with respect to time."""
    return np.gradient(data, time)

def extract_time_range(t, Ip, time_range):
    """Extracts the portion of data within the specified time range."""
    indices = np.where((t >= time_range[0]) & (t <= time_range[1]))
    return t[indices], Ip[indices]

def process_shot_data(shot_group, time_range, window_size):
    """Processes shot data, applies smoothing, and calculates the derivative."""
    Ip = -shot_group['Ip'][:]
    t = shot_group['t'][:]
    t_range, Ip_range = extract_time_range(t, Ip, time_range)
    Ip_smoothed = moving_average(Ip_range, window_size)
    dIp_dt_smoothed = calculate_derivative(Ip_smoothed, t_range)
    return t_range, Ip_smoothed, dIp_dt_smoothed

def plot_all_shots_data(shots_data, ylabel, title, plot_type='Ip'):
    """Plots all shots data on a single figure."""
    plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 16})
    for shot_name, data in shots_data.items():
        if plot_type == 'Ip':
            plt.plot((data['t']-62.525)*1000, data['Ip_smoothed']*1e-6, label=f'{shot_name}', linewidth=2.0)
        elif plot_type == 'dIp/dt':
            plt.plot(data['t']-62.525, data['dIp_dt_smoothed'], label=f'{shot_name}')
    plt.xlabel('$t$ [ms]')
    plt.ylabel(ylabel)
    #plt.title(title)
    plt.legend()
    #plt.grid(True)
    plt.show()

def save_shot_data(shot_name, t, Ip_smoothed, dIp_dt_smoothed, folder='saved_data'):
    """Saves the smoothed Ip and its derivative data to a file for a shot."""
    # Ensure the folder exists
    os.makedirs(folder, exist_ok=True)
    
    # Define file paths
    ip_file_path = os.path.join(folder, f'{shot_name}_Ip_smoothed.txt')
    derivative_file_path = os.path.join(folder, f'{shot_name}_dIp_dt_smoothed.txt')
    
    # Save the data
    np.savetxt(ip_file_path, np.vstack((t, Ip_smoothed)).T, header='Time (s), Smoothed Ip (A)')
    np.savetxt(derivative_file_path, np.vstack((t, dIp_dt_smoothed)).T, header='Time (s), dIp/dt (A/s)')
    print(f'Data saved for {shot_name}.')

def process_and_plot_data(fp, time_range=(62.525, 62.57), window_size=5):
    """Processes, plots, and saves data for all shots."""
    shots_data = {}
    with h5py.File(fp, 'r') as file:
        for shot_name in list(file.keys())[:-1]:
            shot_group = file[shot_name]
            t_range, Ip_smoothed, dIp_dt_smoothed = process_shot_data(shot_group, time_range, window_size)
            print(t_range.size)
            shots_data[shot_name] = {'t': t_range, 'Ip_smoothed': Ip_smoothed, 'dIp_dt_smoothed': dIp_dt_smoothed}
            
            # Save the data for each shot
            save_shot_data(shot_name, t_range, Ip_smoothed, dIp_dt_smoothed)
    
    # Plotting smoothed Ip for all shots
    plot_all_shots_data(shots_data, '$I_{\mathrm{p}}$ [MA]', 'Smoothed Plasma Current as a Function of Time', plot_type='Ip')

    # Plotting derivative of Ip for all shots
    plot_all_shots_data(shots_data, 'dIp/dt (A/s)', 'Derivative of Plasma Current as a Function of Time', plot_type='dIp/dt')

# Example usage
fp = 'JET-IP.h5'  # Adjust this path as necessary
process_and_plot_data(fp)

