import os
import numpy as np
import matplotlib.pyplot as plt
import re

discharge = '85445'
# Base directory where the folders are stored
base_dir = f'../JETresults/parameter_scans/{discharge}'  # Adjust YOUR_DISCHARGE_VALUE accordingly

# Prepare lists to store dBB, assimilation values, and final I_RE values
dBB_list = []
assimilation_list = []
I_RE_final_values = []

# Regular expression to extract dBB and assimilation values from folder names
pattern = re.compile(r"dBB_([0-9.]+)_assimilation_([0-9.]+)")

# Iterate over directories in the base directory
for folder in os.listdir(base_dir):
    match = pattern.match(folder)
    if match:
        dBB, assimilation = map(float, match.groups())
        dBB_list.append(dBB)
        assimilation_list.append(assimilation)

        # Construct the path to the I_RE.txt file
        file_path = os.path.join(base_dir, folder, "I_RE.txt")
        # Load the last value of I_RE data from the file
        with open(file_path, 'r') as file:
            lines = file.readlines()
            I_RE_final = float(lines[-1].strip())  # Read the last line and convert to float
        I_RE_final_values.append(I_RE_final)



# Convert lists to numpy arrays for plotting
dBB_values = np.unique(dBB_list)
assimilation_values = np.unique(assimilation_list)
I_RE_final_matrix = np.zeros((len(dBB_values), len(assimilation_values)))

# Populate the I_RE_final_matrix with I_RE_final_values
for i, dBB in enumerate(dBB_values):
    for j, assimilation in enumerate(assimilation_values):
        index = dBB_list.index(dBB) + assimilation_list.index(assimilation)
        I_RE_final_matrix[i, j] = I_RE_final_values[index]

# Plotting
dBB_mesh, assimilation_mesh = np.meshgrid(assimilation_values, dBB_values)

# Create the contour plot
plt.figure(figsize=(10, 7))
contour = plt.contourf(dBB_mesh, assimilation_mesh, I_RE_final_matrix, levels=20, cmap='plasma')
plt.colorbar(contour, label='Final $I_{re}$ (A)')
plt.title('Final Runaway Electron Current ($I_{re}$)')
plt.xlabel('Assimilation (%)')
plt.ylabel('dBB')
plt.show()
