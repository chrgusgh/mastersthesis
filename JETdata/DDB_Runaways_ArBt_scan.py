import h5py
import numpy as np
import matplotlib.pyplot as plt

Ar_str          = 'Ar_Pam3'
D2_str          = 'D2_Pam3'
Ip_str          = 'Ip_A'
Ne_LIDAR_str    = 'Ne_profile_LIDAR_mm3'
R_Ne_LIDAR_str  = 'R_Ne_profile_LIDAR_m'
R_Te_ECE_str    = 'R_Te_profile_ECE_m'
R_Te_LIDAR_str  = 'R_Te_profile_LIDAR_m'
Te_ECE_str      = 'Te_profile_ECE_eV'
Te_LIDAR_str    = 'Te_profile_LIDAR_eV'

fp = 'DDB_Runaways_ArBt_scan.h5'
f = h5py.File(fp, 'r')

p0 = f['85021']
p1 = f['85445']
p2 = f['85450']
p3 = f['85451']
p4 = f['85453']
p5 = f['85943']

SHOT_nr = 0 # Choose one
SHOT = [p0, p1, p2, p3, p4, p5]

plt.figure(figsize=(10,6))
for SHOT_nr in range(6):

    Ar = SHOT[SHOT_nr][Ar_str][:]
    D2 = SHOT[SHOT_nr][D2_str][:]
    Ip = SHOT[SHOT_nr][Ip_str][:]
    Ne_LIDAR = SHOT[SHOT_nr][Ne_LIDAR_str][:]
    R_Ne_LIDAR = SHOT[SHOT_nr][R_Ne_LIDAR_str][:]
    R_Te_ECE = SHOT[SHOT_nr][R_Te_ECE_str][:]
    R_Te_LIDAR = SHOT[SHOT_nr][R_Te_LIDAR_str][:]
    Te_ECE = SHOT[SHOT_nr][Te_ECE_str][:]
    Te_LIDAR = SHOT[SHOT_nr][Te_LIDAR_str][:]
    
    print('shot:', SHOT_nr)
    print('Ar:', Ar)
    print('D2:', D2)
    #print(R_Ne_LIDAR.shape)
    #print(R_Te_ECE.shape)
    #print(R_Te_LIDAR.shape)

    R0 = 2.96
    LA = R_Te_LIDAR - R0
    #print(LA)
    smt = (LA[-1] - LA[0])/2
    #print(smt)

    plt.plot(R_Te_LIDAR[22:-8] - R0, Te_LIDAR[22:-8], label=f'shot: {SHOT_nr + 1}')

plt.xlabel('r')
plt.ylabel('Te')
plt.legend(loc='best')
plt.title('Te LIDAR')

for SHOT_nr in range(6):

    Ar = SHOT[SHOT_nr][Ar_str][:]
    D2 = SHOT[SHOT_nr][D2_str][:]
    Ip = SHOT[SHOT_nr][Ip_str][:]
    Ne_LIDAR = SHOT[SHOT_nr][Ne_LIDAR_str][:]
    R_Ne_LIDAR = SHOT[SHOT_nr][R_Ne_LIDAR_str][:]
    R_Te_ECE = SHOT[SHOT_nr][R_Te_ECE_str][:]
    R_Te_LIDAR = SHOT[SHOT_nr][R_Te_LIDAR_str][:]
    Te_ECE = SHOT[SHOT_nr][Te_ECE_str][:]
    Te_LIDAR = SHOT[SHOT_nr][Te_LIDAR_str][:]

    R0 = 2.92
    LA = R_Te_ECE - R0
    #print(LA)
    #smt = (LA[-1] - LA[0])/2
    #print(smt)
    plt.figure(figsize=(10,6))
    R_TE_ECE = R_Te_ECE - R0
    print(R_TE_ECE.shape)
    plt.plot(R_TE_ECE[:], Te_ECE, label=f'shot: {SHOT_nr + 1}')
    plt.xlabel('r')
    plt.ylabel('Te')
    plt.title('Te ECE')
    plt.legend(loc='best')

for SHOT_nr in range(6):

    Ar = SHOT[SHOT_nr][Ar_str][:]
    D2 = SHOT[SHOT_nr][D2_str][:]
    Ip = SHOT[SHOT_nr][Ip_str][:]
    Ne_LIDAR = SHOT[SHOT_nr][Ne_LIDAR_str][:]
    R_Ne_LIDAR = SHOT[SHOT_nr][R_Ne_LIDAR_str][:]
    R_Te_ECE = SHOT[SHOT_nr][R_Te_ECE_str][35:-5]
    R_Te_LIDAR = SHOT[SHOT_nr][R_Te_LIDAR_str][:]
    Te_ECE = SHOT[SHOT_nr][Te_ECE_str][35:-5]
    Te_LIDAR = SHOT[SHOT_nr][Te_LIDAR_str][:]

    print(Ne_LIDAR.shape)
    print(R_Ne_LIDAR.shape)
    R0 = 2.92
    plt.figure(figsize=(10,6))
    plt.plot(R_Ne_LIDAR[21:-9] - R0, Ne_LIDAR[21:-9], label=f'shot: {SHOT_nr + 1}')
    plt.xlabel('r')
    plt.ylabel('Ne')
    plt.title('Ne LIDAR')
    plt.legend(loc='best')

plt.show()




