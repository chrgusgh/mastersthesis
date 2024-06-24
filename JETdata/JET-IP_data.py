import numpy as np
import matplotlib.pyplot as plt
import h5py

fp = f'JET-IP.h5'

plt.figure(figsize=(10, 6))

# Open the HDF5 file
with h5py.File(fp, 'r') as file:
    # Iterate over each shot in the file

    keyys = list(file.keys())[:-1]

    for shot_name in keyys:
        print(f"Shot: {shot_name}")
        # Access the group for the current shot
        shot_group = file[shot_name]
        keyz = list(shot_group.keys())
        print(f"Keys in {shot_name}: {keyz}")
        
        # Read the 'Ip' dataset
        Ip = -shot_group['Ip'][:]
        print(f"Ip: {Ip}")
        # Read the 't' dataset
        t = shot_group['t'][:]
        print(f"t: {t}")

        plt.plot(t, Ip, label=f'{shot_name}')

plt.xlabel('Time (t)')
plt.ylabel('Ip (A)')
plt.title('Plasma Current as a Function of Time')
plt.legend()
plt.grid(True)
plt.show()
