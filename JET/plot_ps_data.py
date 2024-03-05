import os
import numpy as np
import matplotlib.pyplot as plt
import re
import utils as utils
from simulation import TokamakSimulation

discharge = '85445'
# Base directory where the folders are stored
base_dir = f'../JETresults/parameter_scans/{discharge}'

# Prepare lists to store dBB, assimilation values, and final I_RE values
dBB_list = []
assimilation_list = []
I_RE_max_values = []
I_Ohm_final_values = []
tau_CQ_values = []

def read_I_RE_final(base_dir, folder):
    # Construct the path to the I_RE.txt file
    file_path_I = os.path.join(base_dir, folder, "I_RE.txt")
    # Load the last value of I_RE data from the file
    with open(file_path_I, 'r') as file:
        return float(file.readlines()[-1].strip())

def read_file(base_dir, folder, filename):
    file_path = os.path.join(base_dir, folder, filename)
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Strip whitespace from the line
            stripped_line = line.strip()
            # Skip lines that start with '#' or are empty
            if stripped_line.startswith('#') or not stripped_line:
                continue
            try:
                # Attempt to convert the line to float and append it to data
                data.append(float(stripped_line))
            except ValueError:
                # If conversion fails, skip the line
                #print(f"Skipping line that cannot be converted to float: {stripped_line}")
                continue
    return data

def read_tau_CQ(base_dir, folder):
    # Path to the tau_CQ.txt file
    file_path_tau_CQ = os.path.join(base_dir, folder, "tau_CQ.txt")

    with open(file_path_tau_CQ, 'r') as file:
        line = file.readline().strip()  # Read and strip the first line
        # Extract the numeric value from the string
        tau_CQ_str = line.split('[')[-1].split(']')[0]  # Assumes the format "tau_CQ: [number]"
        return float(tau_CQ_str)

# Regular expression to extract dBB and assimilation values from folder names
pattern = re.compile(r"dBB_([0-9.]+)_assim_([0-9.]+)")
#simulation = TokamakSimulation()

# Iterate over directories in the base directory
#for folder in os.listdir(base_dir):
#    match = pattern.match(folder)
#    if match:
        #dBB, assimilation = map(float, match.groups())
        #dBB_list.append(dBB)
        #assimilation_list.append(assimilation)

#        I_RE_max = read_file(base_dir, folder, "I_RE.txt")
#        I_RE_max = np.array(I_RE_max)
#        I_RE_max = np.max(I_RE_max)
#        I_RE_max_values.append(I_RE_max)
        
        #I_Ohm = read_file(base_dir, folder, "I_Ohm.txt")
        #I_Ohm = np.array(I_Ohm)

        #t = read_file(base_dir, folder, "t.txt")

        #tau_CQ = red_tau_CQ(base_dir, folder)
        #tau_CQ_values.append(tau_CQ)

        #tau_CQ = utils.calculate_t_CQ(I_Ohm, simulation.Ip0, t)

        #print('dBB:', dBB, 'assim:', assimilation, 'tau_CQ:', tau_CQ)

dBBs = []
assimilations = []

for folder in os.listdir(base_dir):
    match = pattern.match(folder)
    if match:
        dBB, assimilation = match.group(1), match.group(2)
        dBBs.append(dBB)
        assimilations.append(assimilation)

dBBs = np.unique(dBBs)
assimilations = np.unique(assimilations)
I_RE_matrix = np.zeros((len(dBBs), len(assimilations)))
for i in range(len(dBBs)):
    for j in range(len(assimilations)):
        dir_name = f'dBB_{dBBs[i]}_assim_{assimilations[j]}'
        file_path = os.path.join(base_dir, dir_name, "I_RE.txt")
        data = []
        with open(file_path, 'r') as file:
            for line in file:
                # Strip whitespace from the line
                stripped_line = line.strip()
                # Skip lines that start with '#' or are empty
                if stripped_line.startswith('#') or not stripped_line:
                    continue
                try:
                    # Attempt to convert the line to float and append it to data
                    data.append(float(stripped_line))
                except ValueError:
                    # If conversion fails, skip the line
                    #print(f"Skipping line that cannot be converted to float: {stripped_line}")
                    continue
                
        I_RE_matrix[i, j] = np.max(np.array(data))

dBBs_arr = dBBs.astype(float)
assimilations_arr = assimilations.astype(float)
print(I_RE_matrix.shape)
D, A = np.meshgrid(dBBs_arr, assimilations_arr)

# Convert lists to numpy arrays for plotting
#dBB_values = np.array(dBB_list)
#assimilation_values = np.array(assimilation_list)
#I_RE_max_matrix = np.zeros((len(dBB_values), len(assimilation_values)))

# Populate the I_RE_final_matrix with I_RE_final_values
#for i, dBB in enumerate(dBB_values):
#    for j, assimilation in enumerate(assimilation_values):
#        index = dBB_list.index(dBB) + assimilation_list.index(assimilation)
#        I_RE_max_matrix[i, j] = I_RE_max_values[index]

#print(dBB_values, assimilation_values, I_RE_final_matrix.shape)

# Plotting
#dBB_mesh, assimilation_mesh = np.meshgrid(assimilation_values), dBB_values))

# Create the contour plot
plt.figure(figsize=(10, 7))
contour = plt.contourf(D, A, I_RE_matrix.T, cmap='plasma')
plt.colorbar(contour, label='$I_{re}$ (A)')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')

plt.title('Runaway Electron Current ($I_{re}$)')
plt.ylabel('Assimilation')
plt.xlabel('dBB')
plt.show()
