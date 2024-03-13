import os
import numpy as np
import matplotlib.pyplot as plt
import re

def get_dBB_and_assimilation_from_folder_name(folder, pattern):
    match = pattern.match(folder)
    if match:
        return match.group(1), match.group(2)
    return None, None

def read_data_from_file(file_path):
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
    D, A = np.meshgrid(assimilation_values, dBB_values)
    plt.figure(figsize=(10, 7))
    contour = plt.contourf(D, A, I_RE_matrix.T, cmap='plasma')
    plt.colorbar(contour, label='$I_{RE}$ (A)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Runaway Electron Current ($I_{RE}$)')
    plt.xlabel('Assimilation')
    plt.ylabel('dBB')
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

if __name__ == '__main__':
    main()

