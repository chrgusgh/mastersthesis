import numpy as np
import scipy.constants as const

def impurity_density(impurity, Ti, plasma_volume):
    """
    Computes impurity density from the gas law

    Parameters:
    - impurity: Injected impurity amount in units Pa m^3
    - Ti: Initial impurity temperature K
    - plasma_volume: Plasma volume in m^3
    
    Returns:
    - n_i: Impurity density for ion species i.
    """
    N = impurity / (const.k * Ti)
    n_i = N / plasma_volume
    return n_i

def calculate_v_th(T_final, m_e_eV):
    """
    Calculate the thermal velocity of electrons at T = T_final (100 eV)
    
    Parameters:
    - T_final: Electron temperature in eV
    - m_e_eV: Electron mass in eV

    Returns:
    - v_th: Thermal electron velocity in m/s
    """
    return np.sqrt(2 * T_final / m_e_eV)

def calculate_dBB(a, R0, T_final, t_TQ, v_th):
    """
    Calculate the magnetic perturbation parameter dBB.

    Parameters:
    - a: Minor radius of the tokamak.
    - R0: Major radius of the tokamak.
    - t_TQ: Thermal quench time in seconds.
    - v_th: Thermal velocity.

    Returns:
    - dBB: Magnetic perturbation parameter.
    """
    return np.sqrt(a / np.sqrt(const.pi * t_TQ * R0 * v_th))

def calculate_D(a, t_TQ):
    """
    Calculate the diffusion coefficient.

    Parameters:
    - a: Minor radius of tokamak
    - t_TQ: Thermal quench time in seconds.

    Returns:
    - D = Drr: Diffusion coefficient
    """
    return a ** 2 / t_TQ

def calculate_Drr(R0, q, dBB_grid):
    """
    Calculate the diffusion coefficient D based on the Rechester-Rosenbluth formula.

    Parameters:
    - R0: Major radius of the tokamak.
    - q: Charge quantity.
    - c: Speed of light in vacuum.
    - dBB_grid: Grid of magnetic perturbation values.

    Returns:
    - D: Diffusion coefficient.
    """
    return const.pi * R0 * q * const.c * dBB_grid ** 2

def calculate_tau_TQ(t_TQ, T_init, T_final):
    """
    Calculates the thermal quench time variable tau_TQ, from the exponential decay formula.

    Parameters:
    - t_TQ: Thermal quench time in seconds
    - T_init: Initial plasma temperature in eV
    - T_final: Plasma temperature at the end of the TQ, set to 100 eV

    Returns:
    - tau_TQ: Thermal quench time variable in seconds.
    """
    return t_TQ * np.log(T_init / T_final)


