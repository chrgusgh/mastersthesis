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

def calculate_dBB(a, R0, t_TQ, v_th):
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
    return a / np.sqrt(const.pi * t_TQ * R0 * v_th)

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

def calculate_Drr(R0, q, dBB):
    """
    Calculate the diffusion coefficient D based on the Rechester-Rosenbluth formula.

    Parameters:
    - R0: Major radius of the tokamak.
    - q: Charge quantity.
    - c: Speed of light in vacuum.
    - dBB: Magnetic perturbation values.

    Returns:
    - D: Diffusion coefficient.
    """
    return const.pi * R0 * q * const.c * dBB ** 2

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

def calculate_t_CQ(I_Ohm, Ip0, t):
    """
    Calculates the current quench time from the formula.

    Parameters:
    - I_Ohm: Ohmic current array
    - Ip0: Initial plasma current
    - t: Time array for the current quench phase

    Returns:
    - t_CQ: Current quench time in seconds. If calculation fails, returns 0.
    """
    try:
        I_Ohm_init = I_Ohm[0]
        lb = 0.799 * I_Ohm_init
        ub = 0.801 * I_Ohm_init
        indices = np.where((I_Ohm >= lb) & (I_Ohm <= ub))[0]
        if indices.size == 0:
            raise IndexError("No valid index found.")
        I_Ohm_80_idx = indices[0]
        t_80 = t[I_Ohm_80_idx]
        t_last = t[-1]
        # Note minus sign in front
        return -(t_80 - t_last) / (0.8 - I_Ohm[-1] / Ip0)
    except IndexError as e:
        print(f"Error calculating t_CQ: {e}")
        return 0
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 0

