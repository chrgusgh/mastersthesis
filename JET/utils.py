import numpy as np
import scipy.constants as const

def calculate_dBB(a, R0, T_final, t_TQ, v_th):
    """
    Calculate the magnetic perturbation parameter dBB.

    Parameters:
    - a: Minor radius of the tokamak.
    - R0: Major radius of the tokamak.
    - T_final: Final temperature in eV.
    - t_TQ: Thermal quench time in seconds.
    - v_th: Thermal velocity.

    Returns:
    - dBB: Magnetic perturbation parameter.
    """
    return np.sqrt(a / np.sqrt(const.pi * t_TQ * R0 * v_th))

def calculate_D(R0, q, dBB_grid):
    """
    Calculate the diffusion coefficient D based on the magnetic perturbation.

    Parameters:
    - R0: Major radius of the tokamak.
    - q: Charge quantity.
    - c: Speed of light in vacuum.
    - dBB_grid: Grid of magnetic perturbation values.

    Returns:
    - D: Diffusion coefficient.
    """
    return const.pi * R0 * q * const.c * dBB_grid ** 2

def set_transport_conditions(ds, D_st_grid, dBB_st_grid, r_grid, tau_TQ_grid, MODE):
    """
    Set transport conditions for the DREAM settings object based on the specified parameters.

    Parameters:
    - ds: DREAM settings object to configure.
    - D_st_grid: 2D grid of diffusion coefficients.
    - dBB_st_grid: 2D grid of magnetic perturbation parameters.
    - r_grid: Radial grid.
    - tau_TQ_grid: Time grid for thermal quench.
    - MODE: Operational mode for the simulation.
    """
    ds.eqsys.n_re.transport.setBoundaryCondition(driver.Transport.BC_F_0)
    ds.eqsys.n_re.transport.prescribeDiffusion(D_st_grid, r=r_grid, t=tau_TQ_grid)

    ds.eqsys.T_cold.transport.setMagneticPerturbation(dBB=dBB_st_grid, r=r_grid, t=tau_TQ_grid)
    ds.eqsys.T_cold.transport.setBoundaryCondition(driver.Transport.BC_F_0)

    if MODE in (driver.MODE_ISOTROPIC, driver.MODE_KINETIC):
        ds.eqsys.f_hot.transport.setBoundaryCondition(driver.Transport.BC_F_0)
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB=dBB_st_grid, r=r_grid, t=tau_TQ_grid)


def set_magnetic_perturbation(ds, MODE, nr, nt, T0, a, R0):
    """
    Sets the simulation settings for the DREAM object based on the specified parameters.

    Parameters:
    - ds: DREAM settings object to configure.
    - MODE: Operational mode for the simulation.
    - nr: Number of radial points.
    - nt: Number of time points.
    """
    q = 1
    m_e_eV = const.m_e * const.c**2 / const.e
    T_init = T0[0]  # eV, initial temperature
    T_final = 100  # eV, final temperature
    t_TQ = 0.0002  # s, thermal quench time
    tau_TQ = t_TQ * np.log(T_init / T_final)
    tau_TQ_grid = np.linspace(0, tau_TQ, nt)
    r_grid = np.linspace(0, a, nr)
    v_th = np.sqrt(2 * T_final / m_e_eV)  # Thermal velocity at T_final
    
    dBB = calculate_dBB(a, R0, T_final, t_TQ, v_th)
    dBB_grid = np.ones(nr) * dBB
    ones_2D = np.ones((nt, nr))
    dBB_st_grid = dBB_grid * ones_2D
    
    D = calculate_D(R0, q, const.c, dBB_grid)
    D_st_grid = D * ones_2D

    set_transport_conditions(ds, D_st_grid, dBB_st_grid, r_grid, tau_TQ_grid, MODE)

