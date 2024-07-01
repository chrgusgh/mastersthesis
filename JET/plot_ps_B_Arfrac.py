import os
import numpy as np
import matplotlib.pyplot as plt
import re

def get_B_factor_and_Arfrac_from_folder_name(folder, pattern):
    """
    Extract B_factor and Arfrac values from the folder name using regex.
    
    Parameters:
    - folder (str): Folder name containing the B_factor and Arfrac values.
    - pattern (re.Pattern): Compiled regex pattern to match the folder name format.
    
    Returns:
    - (str, str): Tuple containing B_factor and Arfrac values as strings. Returns (None, None) if no match is found.
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

def create_meshgrid_for_contour_plot(B_factors, Arfracs, base_dir, pattern):
    """
    Create a meshgrid and populate I_RE_matrix for contour plotting based on B_factor and Arfrac.
    
    Parameters:
    - B_factors (list): List of B_factor values as strings.
    - Arfracs (list): List of Arfrac values as strings.
    - base_dir (str): Base directory path containing folders for each parameter combination.
    - pattern (re.Pattern): Compiled regex pattern to match the folder names.
    
    Returns:
    - (np.ndarray, np.ndarray, np.ndarray): Tuple containing arrays for B_factor values, Arfrac values, and the I_RE matrix.
    """
    B_factor_values = np.unique(np.array(B_factors).astype(float))
    Arfrac_values = np.unique(np.array(Arfracs).astype(float))
    I_RE_matrix = np.zeros((len(B_factor_values), len(Arfrac_values)))
    
    for folder in os.listdir(base_dir):
        B_factor, Arfrac = get_B_factor_and_Arfrac_from_folder_name(folder, pattern)
        if B_factor is not None and Arfrac is not None:
            B_factor = float(B_factor)
            Arfrac = float(Arfrac)
            i = np.where(B_factor_values == B_factor)[0][0]
            j = np.where(Arfrac_values == Arfrac)[0][0]
            file_path = os.path.join(base_dir, folder, "I_RE.txt")
            I_RE_data = read_data_from_file(file_path)
            I_RE_matrix[i, j] = np.max(I_RE_data)
            
    return B_factor_values, Arfrac_values, I_RE_matrix

def plot_contour(B_factor_values, Arfrac_values, I_RE_matrix):
    """
    Plot a contour map showing the relationship between B_factor, Arfrac, and I_RE.
    
    Parameters:
    - B_factor_values (np.ndarray): Array of B_factor values.
    - Arfrac_values (np.ndarray): Array of Arfrac values.
    - I_RE_matrix (np.ndarray): 2D array containing I_RE values corresponding to each B_factor and Arfrac combination.
    """
    A, D = np.meshgrid(np.sort(Arfrac_values), np.sort(B_factor_values))
    plt.figure(figsize=(10, 7))
    plt.rcParams.update({'font.size': 20})
    #levels = np.linspace(I_RE_matrix.min(), I_RE_matrix.max(), 50)  # Adjust 50 to increase/decrease the number of levels
    contour = plt.contourf(A, D, I_RE_matrix*1e-6, cmap='plasma', levels=28)
    #contour = plt.contourf(A, D, I_RE_matrix*1e-6, cmap='plasma', levels= 100)
    plt.colorbar(contour, label='$I_{\mathrm{RE,\, max}}$ [MA]')
    for collection in contour.collections:
        collection.set_rasterized(True)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.title('Maximum Runaway Electron Current ($I_{re, max}$)')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('Ar fraction')
    plt.ylabel('$C_B$')
    plt.show()

def main():
    discharge = '85943'
    base_dir = '../../../../../mnt/DISK4/christiang/resultat/parameter_scans/B_factor_Arfrac_scans/85943_the_true_cD_3.5'
    #base_dir = '../../../../../mnt/DISK4/christiang/resultat/parameter_scans/B_factor_Arfrac_scans/85943_the_one'
    pattern = re.compile(r"B_factor_([0-9.]+)_Arfrac_([0-9.]+)")
    
    B_factors, Arfracs = [], []
    for folder in os.listdir(base_dir):
        B_factor, Arfrac = get_B_factor_and_Arfrac_from_folder_name(folder, pattern)
        if B_factor and Arfrac:
            B_factors.append(B_factor)
            Arfracs.append(Arfrac)
    
    B_factor_values, Arfrac_values, I_RE_matrix = create_meshgrid_for_contour_plot(B_factors, Arfracs, base_dir, pattern)
    
    plot_contour(B_factor_values, Arfrac_values, I_RE_matrix)

if __name__ == '__main__':
    main() 
