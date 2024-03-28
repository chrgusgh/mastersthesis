import os
import numpy as np
import matplotlib.pyplot as plt
import re

def get_dBB_and_assimilation_from_folder_name(folder, pattern):
    """
    Extract dBB and assimilation values from the folder name using regex.
    
    Parameters:
    - folder (str): Folder name containing the dBB and assimilation values.
    - pattern (re.Pattern): Compiled regex pattern to match the folder name format.
    
    Returns:
    - (str, str): Tuple containing dBB and assimilation values as strings. Returns (None, None) if no match is found.
    """
    match = pattern.match(folder)
    if match:
        return match.group(1), match.group(2)
    return None, None

def read_data_from_file(file_path):
    """
    Read and return numerical data from a file, ignoring comments and empty lines.
    
    Parameters:
    - file_path (str): Path to the file containing the data.
    
    Returns:
    - np.ndarray: Array of data read from the file.
    """
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()
            if stripped_line.startswith('#') or not stripped_line:
                continue
            try:
                data.append(float(stripped_line))
            except ValueError:
                continue
    return np.array(data)

def create_meshgrid_for_contour_plot(dBBs, assimilations, base_dir, pattern):
    """
    Create a meshgrid and populate I_RE_matrix for contour plotting based on dBB and assimilation.
    
    Parameters:
    - dBBs (list): List of dBB values as strings.
    - assimilations (list): List of assimilation values as strings.
    - base_dir (str): Base directory path containing folders for each parameter combination.
    - pattern (re.Pattern): Compiled regex pattern to match the folder names.
    
    Returns:
    - (np.ndarray, np.ndarray, np.ndarray): Tuple containing arrays for dBB values, assimilation values, and the I_RE matrix.
    """
    dBB_values = np.unique(np.array(dBBs).astype(float))
    assimilation_values = np.unique(np.array(assimilations).astype(float))
    I_RE_matrix = np.zeros((len(dBB_values), len(assimilation_values)))
    
    for folder in os.listdir(base_dir):
        dBB, assimilation = get_dBB_and_assimilation_from_folder_name(folder, pattern)
        if dBB is not None and assimilation is not None:
            dBB = float(dBB)
            assimilation = float(assimilation)
            i = np.where(dBB_values == dBB)[0][0]
            j = np.where(assimilation_values == assimilation)[0][0]
            file_path = os.path.join(base_dir, folder, "I_RE.txt")
            I_RE_data = read_data_from_file(file_path)
            I_RE_matrix[i, j] = np.max(I_RE_data)
            
    return dBB_values, assimilation_values, I_RE_matrix

def plot_contour(dBB_values, assimilation_values, I_RE_matrix):
    """
    Plot a contour map showing the relationship between dBB, assimilation, and I_RE.
    
    Parameters:
    - dBB_values (np.ndarray): Array of dBB values.
    - assimilation_values (np.ndarray): Array of assimilation values.
    - I_RE_matrix (np.ndarray): 2D array containing I_RE values corresponding to each dBB and assimilation combination.
    """
    A, D = np.meshgrid(assimilation_values, dBB_values)
    #print(A)
    #print(D)
    plt.figure(figsize=(10, 7))
    contour = plt.contourf(A, D, I_RE_matrix, cmap='plasma')
    plt.colorbar(contour, label='$I_{RE}$ (A)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Runaway Electron Current ($I_{RE}$)')
    plt.xlabel('Assimilation')
    plt.ylabel('$\delta B / B$')
    plt.show()

def main():
    discharge = '85445'
    base_dir = f'../JETresults/parameter_scans/{discharge}'
    pattern = re.compile(r"dBB_([0-9.]+)_assim_([0-9.]+)")
    
    dBBs, assimilations = [], []
    for folder in os.listdir(base_dir):
        dBB, assimilation = get_dBB_and_assimilation_from_folder_name(folder, pattern)
        if dBB and assimilation:
            dBBs.append(dBB)
            assimilations.append(assimilation)
    
    dBB_values, assimilation_values, I_RE_matrix = create_meshgrid_for_contour_plot(dBBs, assimilations, base_dir, pattern)
    plot_contour(dBB_values, assimilation_values, I_RE_matrix)
    #print(I_RE_matrix.shape)
    #print(I_RE_matrix)
if __name__ == '__main__':
    main()

