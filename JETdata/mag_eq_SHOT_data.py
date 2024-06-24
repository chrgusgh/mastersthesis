import h5py
import numpy as np
import matplotlib.pyplot as plt

SHOT = 0
SHOTS = ['85021','85445','85450','85451','85453','85943']
fp = f'mag_eq_{SHOTS[SHOT]}_data.h5'
equil = h5py.File(fp, 'r')['equil']

eq_list = list(equil)
#print(eq_list)

ent = 'ptx' 
#print(ent, ':')
entry = equil[ent][:]
#print(entry[:])

radii_indices = [2]
ptBx = equil['ptBx'][:]
ptBy = equil['ptBy'][:]
ptBPhi = equil['ptBPHI'][:]
theta = equil['theta'][:]

print(ptBx)
print(ptBy)
print(ptBPhi)
# Calculate B for each radius and theta
B = np.sqrt(ptBx**2 + ptBy**2 + ptBPhi**2)

# Plot B as a function of theta for the selected radii
plt.figure(figsize=(10, 6))

for idx in range(80):
    B_idx = B[:, idx]
    #plt.plot(theta, B_idx, label=f'Radius index {idx}')   
    #plt.xlabel('Theta (radians)')
    #plt.ylabel('B')
    #plt.title('Magnetic Field Strength (B) as a Function of Theta for Selected Radii')
    #plt.legend()
    #plt.show()

    #print(idx)

    unique_elements, counts = np.unique(B_idx, return_counts=True)
    ue_duplicates = unique_elements[counts > 1]
    #print(ue_duplicates)

    smallest_element = np.min(B_idx)
    largest_element = np.max(B_idx)
    
    tol = 1e-10

    #if len(ue_duplicates) > 0:
        #if np.abs(ue_duplicates[0] - smallest_element) < tol or np.abs(ue_duplicates[0] - largest_element) < tol:
            #print('radius index:', idx)
            #print('thetas:', theta)
    #print('B_idx:', B_idx)
            #print(B_idx[-1]- B_idx[0])
            #print('duplicates B:', duplicates)
            #print('smallest B:', smallest_element)
            #print('largest B:', largest_element)
