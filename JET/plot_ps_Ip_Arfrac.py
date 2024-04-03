import os
import numpy as np
import matplotlib.pyplot as plt
import re

def get_Ip_factor_and_Arfrac_from_folder_name(folder, pattern):
    """
    Extract Ip_factor and Arfrac values from the folder name using regex.
    
    Parameters:
    - folder (str): Folder name containing the Ip_factor and Arfrac values.
    - pattern (re.Pattern): Compiled regex pattern to match the folder name format.
    
    Returns:
    - (str, str): Tuple containing Ip_factor and Arfrac values as strings. Returns (None, None) if no match is found.
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

def create_meshgrid_for_contour_plot(Ip_factors, Arfracs, base_dir, pattern):
    """
    Create a meshgrid and populate I_RE_matrix for contour plotting based on Ip_factor and Arfrac.
    
    Parameters:
    - Ip_factors (list): List of Ip_factor values as strings.
    - Arfracs (list): List of Arfrac values as strings.
    - base_dir (str): Base directory path containing folders for each parameter combination.
    - pattern (re.Pattern): Compiled regex pattern to match the folder names.
    
    Returns:
    - (np.ndarray, np.ndarray, np.ndarray): Tuple containing arrays for Ip_factor values, Arfrac values, and the I_RE matrix.
    """
    Ip_factor_values = np.unique(np.array(Ip_factors).astype(float))
    Arfrac_values = np.unique(np.array(Arfracs).astype(float))
    I_RE_matrix = np.zeros((len(Ip_factor_values), len(Arfrac_values)))
    
    counter = 0
    for folder in os.listdir(base_dir):
        Ip_factor, Arfrac = get_Ip_factor_and_Arfrac_from_folder_name(folder, pattern)
        if Ip_factor is not None and Arfrac is not None:
            counter = counter + 1
            Ip_factor = float(Ip_factor)
            Arfrac = float(Arfrac)
            #print(Ip_factor)
            #print(Arfrac)
            #print(counter)
            i = np.where(Ip_factor_values == Ip_factor)[0][0]
            j = np.where(Arfrac_values == Arfrac)[0][0]
            #print('i:',i)
            #print('j:', j)
            file_path = os.path.join(base_dir, folder, "I_RE.txt")
            I_RE_data = read_data_from_file(file_path)
            I_RE_matrix[i, j] = np.max(I_RE_data)
            
    return Ip_factor_values, Arfrac_values, I_RE_matrix

def plot_contour(Ip_factor_values, Arfrac_values, I_RE_matrix):
    """
    Plot a contour map showing the relationship between Ip_factor, Arfrac, and I_RE.
    
    Parameters:
    - Ip_factor_values (np.ndarray): Array of Ip_factor values.
    - Arfrac_values (np.ndarray): Array of Arfrac values.
    - I_RE_matrix (np.ndarray): 2D array containing I_RE values corresponding to each Ip_factor and Arfrac combination.
    """
    A, D = np.meshgrid(Arfrac_values, Ip_factor_values)
    plt.figure(figsize=(10, 7))
    levels = np.linspace(I_RE_matrix.min(), I_RE_matrix.max(), 50)  # Adjust 50 to increase/decrease the number of levels
    contour = plt.contourf(A, D, I_RE_matrix, cmap='plasma')
    plt.colorbar(contour, label='$I_{re, max}$ (A)')
    plt.title('Maximum Runaway Electron Current ($I_{re, max}$)')
    plt.xlabel('Arfrac')
    plt.ylabel('Ip_factor')
    plt.show()

def main():
    discharge = '85943'
    base_dir = f'../../../../../mnt/DISK4/christiang/resultat/parameter_scans/Ip_Arfrac_scans/{discharge}_the_one'
    pattern = re.compile(r"Ip_factor_([0-9.]+)_Arfrac_([0-9.]+)")
    
    Ip_factors, Arfracs = [], []
    for folder in os.listdir(base_dir):
        Ip_factor, Arfrac = get_Ip_factor_and_Arfrac_from_folder_name(folder, pattern)
        if Ip_factor and Arfrac:
            Ip_factors.append(Ip_factor)
            Arfracs.append(Arfrac)
    
    Ip_factor_values, Arfrac_values, I_RE_matrix = create_meshgrid_for_contour_plot(Ip_factors, Arfracs, base_dir, pattern)
    plot_contour(Ip_factor_values, Arfrac_values, I_RE_matrix)

if __name__ == '__main__':
    main()

