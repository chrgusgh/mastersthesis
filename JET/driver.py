# General routines for running a TCV disruption simulation

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from pathlib import Path
from scipy.constants import c, e
from scipy.constants import k as k_B
import sys
import utils as utils
import inputs as inputs

#sys.path.append('/mnt/HDD/software/DREAM')
#sys.path.append('/mnt/HDD/software/DREAM/build/dreampyface/cxx')
sys.path.append('../../DREAM/py/')
from DREAM import DREAMSettings, runiface
#import dreampyface

import DREAM.Settings.AdvectionInterpolation as AI
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as OhmicCurrent
import DREAM.Settings.Equations.RunawayElectrons as RunawayElectrons
import DREAM.Settings.Equations.RunawayElectronDistribution as FRe
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TimeStepper as TimeStepper
import DREAM.Settings.TransportSettings as Transport


# TODO: Update
R0 = inputs.R0
a  = inputs.a
#b  = 0.30
#B0 = 1.45

MODE_FLUID = 1
MODE_ISOTROPIC = 2
MODE_KINETIC = 3

SCAN_NONE = 0
SCAN_NEON = 1
SCAN_NRE  = 2
SCAN_CURRENT = 3
SCAN_TEMPERATURE = 4
SCAN_TEMPERATURE_PROFILE = 5
SCAN_DRR = 6
SCAN_D2 = 7
SCAN_VLOOP = 8


def modename(mode):
    if mode == MODE_FLUID: return 'fluid'
    elif mode == MODE_ISOTROPIC: return 'isotropic'
    elif mode == MODE_KINETIC: return 'kinetic'
    else: return ''


def longsetname(prefix, mode, phase, io='settings', extension=''):
    """
    Generate a settings file name.
    """
    return longoutname(prefix=prefix, mode=mode, phase=phase, io=io, extension=extension)


def longoutname(prefix, mode, phase, io='output', extension=''):
    """
    Generate an output file name.
    """
    base = f'{prefix}_{modename(mode)}_{io}_{phase}'

    filename = f'{base}{extension}.h5'
    p = Path(filename).parent.resolve()

    if not p.exists():
        p.mkdir(parents=True)

    return filename


def bar(pressure, frac_assim=0.5, T=300):
    """
    Convert the given pressure (in bar) to a particle density. The
    pressure is related to the gas temperature T and density n via

      p = n*T
    
    Since 1 bar = 100 kPa we have that the atomic density at a pressure
    p (given in bar) is

                       100 000
      (n [m^-3]) = -------------- * (p [bar])
                    k_B * (T [J])

    :param pressure: Pressure in bar to convert.
    :param frac_assim: Fraction of gas that is assimilated in the plasma.
    :param T: Gas temperature (kelvin).
    """
    #return 1e5 * pressure / (k_B * T)
    # TODO implement actual conversion
    return 1e-2 * pressure / (k_B*T)


def get_profile(name, x, y, default, evaluate=False):
    """
    Ensures that the given input profile y(x) is valid.

    :param x:        Radial grid for the profile.
    :param y:        Plasma profile.
    :param evaluate: If ``True``, returns a numpy array rather than a callable.
    """

    if x is None and y is None:
        x = np.linspace(0, a)
        y = default
    elif x is None and callable(y):
        x = np.linspace(0, a)
    elif x is None and np.isscalar(y):
        if type(y) == float or type(y) == np.float64:
            y = np.array([y])
        x = 0
    elif len(x) != len(y):
        raise Exception(f"Invalid specification of {name} profile.")
    elif not callable(y) and not evaluate:
        y = UnivariateSpline(x, y)
    if evaluate and callable(y):
        return x, y(x)
    else:
        return x, y


def terminate_current(Ip, sim):
    """
    Termination function which stops time stepping when the plasma
    current has reached the target value.
    """
    Ip_data = sim.unknowns.getData('I_p')

    # Return True when simulation should terminate
    return Ip_data['x'][-1,0] >= Ip


def setup_runawaygrid(ds, equilibrium):
    """
    Set up the runaway grid for the given simulation object.
    """
    ds.runawaygrid.setEnabled(True)
    #ds.runawaygrid.setPmax(25)
    ds.runawaygrid.setPmax(50)
    #ds.runawaygrid.setPmin(PMAX)
    ds.runawaygrid.setNp(100)
    ds.eqsys.f_re.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_TCDF, ad_jac=DistFunc.AD_INTERP_JACOBIAN_UPWIND)
    ds.runawaygrid.setBiuniformGrid(psep=5, npsep=50)

    # Toroidal or cylindrical geometry?
    if equilibrium:
        ds.runawaygrid.setTrappedPassingBoundaryLayerGrid(dxiMax=0.1)
    else:
        ds.runawaygrid.setNxi(30)
        #ds.runawaygrid.setBiuniformGrid(thetasep=0.6, nthetasep_frac=0.5)
        ds.runawaygrid.setBiuniformGrid(thetasep=np.pi-0.6, nthetasep_frac=0.5)


def generate_baseline(mode=MODE_ISOTROPIC, equilibrium=None, nt=1, tMax=1e-11, nr=10, reltol=1e-6, nre0=None, nre0_r=0, j0=None, j0r=None, T0=1e3, T0r=None, tauwall=None, Vloop=None, withfre=False, E0=None, E0r=None, Ip0=200e3, n0=1e19, n0r=None, nxi_hot=15, np1_hot=80, np2_hot=60, pmax_hot=0.8, dBB0 = 1e-3, runInit=True, verboseInit=False, prefix='output/generic', extension='', **kwargs):
    """
    Generate a baseline TCV disruption simulation object.

    :param mode:        Type of hot electron model to use (fluid or isotropic).
    :param equilibrium: Name of equilibrium data file to use, or ``None`` for generic TCV-like cylindrical equilibrium.
    """
    j0r, j0 = get_profile('current density', j0r, j0, default=lambda r : (1 - (1-0.001**(1/0.41))*(r/a)**2)**0.41)
    T0r, T0 = get_profile('electron temperature', T0r, T0, default=lambda x : 1e3, evaluate=True)

    setname = lambda phase : longsetname(mode=mode, phase=phase, prefix=prefix, extension=extension)
    outname = lambda phase : longoutname(mode=mode, phase=phase, prefix=prefix, extension=extension)

    ds = DREAMSettings()

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONLESS

    #############################################
    # STEP 1 -- Calculate initial electric field
    #############################################

    if equilibrium:
        ds.radialgrid.setNumerical(equilibrium)
    else:
        # Cylindrical TCV-like equilibrium
        ds.radialgrid.setType(RadialGrid.TYPE_CYLINDRICAL)
        #ds.radialgrid.setMajorRadius(R0)
        ds.radialgrid.setMinorRadius(a)
        ds.radialgrid.setB0(B0)

    ds.radialgrid.setWallRadius(1.1*a) # An assumption
    ds.radialgrid.setNr(nr)

    ds.eqsys.E_field.setType(ElectricField.TYPE_PRESCRIBED_OHMIC_CURRENT)
    ds.eqsys.j_ohm.setCurrentProfile(j0(j0r), radius=j0r, Ip0=Ip0)
    # Use Sauter conductivity
    ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)
    
    # Disable kinetic grids during conductivity simulation
    ds.runawaygrid.setEnabled(False)
    ds.hottailgrid.setEnabled(False)

    # Prescribe initial runaway seed?
    if (nre0 is not None) and (nre0 != 0):
        ds.eqsys.n_re.setInitialProfile(nre0, radius=nre0_r)

        if withfre:
            setup_runawaygrid(ds, equilibrium)
            ds.eqsys.f_re.setInitType(FRe.INIT_XI_NEGATIVE)

    # The initial temperature profile only enters in the shape
    # of the initial Maxwellian on the hot grid. This temperature
    # is that of the cold population in the TQ simulation.
    #ds.eqsys.T_cold.setPrescribedData(T0, radius=T0r)
    if mode == MODE_FLUID or mode == MODE_KINETIC:
        ds.eqsys.T_cold.setPrescribedData(temperature=T0, radius=T0r)
    else:
        ds.eqsys.T_cold.setPrescribedData(temperature=1)

    # Set main ion density
    #if not np.isscalar(T0r):
    #    n0 *= np.ones(T0r.shape)

    n0r, n0 = get_profile('electron density', n0r, n0, default=1e19)
    
    ds.eqsys.n_i.addIon(name='D', Z=1, Z0=1, iontype=Ions.IONS_DYNAMIC, T=T0, r=T0r, n=n0(T0r))

    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(nt)

    ds.other.include('fluid', 'scalar')

    do = runiface(ds, 'output.h5', quiet=False)

    ###############################################
    # STEP 2 -- Prescribe initial current density
    ###############################################
    nfree, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    if mode != MODE_FLUID:
        ds.hottailgrid.setEnabled(True)
        ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

        if mode == MODE_ISOTROPIC:
            ds.hottailgrid.setNxi(1)
        elif equilibrium:
            ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=2/nr)
        else:
            ds.hottailgrid.setNxi(nxi_hot)

        ds.hottailgrid.setPmax(pmax_hot)
        if mode == MODE_KINETIC:
            pth = np.sqrt(2*20/511e3)
            ds.hottailgrid.setNp(np1_hot+np2_hot)
            ds.hottailgrid.setBiuniformGrid(psep=7*pth, npsep=np1_hot)
            ds.eqsys.j_ohm.setCorrectedConductivity(True)
        else:
            ds.hottailgrid.setNp(np1_hot)

        ds.eqsys.j_ohm.setCorrectedConductivity(False)

        nre = 0
        if nre0 is not None:
            nre = nre0

        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=((1-1e-3)*nfree - nre)*np.ones(rn0.shape), T0=T0, rT0=T0r)

        #ds.eqsys.f_hot.setHotRegionThreshold(5)

        ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=AI.AD_INTERP_TCDF, ad_jac=AI.AD_INTERP_JACOBIAN_UPWIND)
        # Do not include Jacobian elements for d f_hot / d n_i, i.e.
        # derivatives with respect to the ion densities (which would take
        # up *significant* space in the matrix)
        ds.eqsys.f_hot.enableIonJacobian(False)

    # Set electric field profile
    # (we set E instead of j here to allow also f_hot to be properly initialized)
    #ds.eqsys.E_field.setPrescribedData(do.eqsys.E_field[-1,:], radius=do.grid.r)
    Einit = do.eqsys.E_field[-1,:]
    Einit_r = do.grid.r[:]

    do.close()

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
    ds.solver.setMaxIterations(100)
    ds.solver.tolerance.set(reltol=reltol)

    if mode != MODE_FLUID:
        ds.solver.tolerance.set('n_hot', abstol=1e5)
        ds.solver.tolerance.set('j_hot', abstol=1)

    # Generate current profile
    INITFILE = outname('init')
    if mode == MODE_KINETIC and runInit:   # No need to run this bit when the bulk is a fluid...
        print('Current simulation')

        # Steady state ohmic current (let t -> inf)
        ds.timestep.setTmax(1)
        ds.timestep.setNt(3)

        ds.eqsys.E_field.setPrescribedData(Einit, radius=Einit_r)

        ds.solver.setVerbose(verboseInit)
        #ds.solver.setDebug(savesystem=True, timestep=1, iteration=0)
        ds.save(setname('init'))
        runiface(ds, INITFILE, quiet=(not verboseInit))

    ###############################################
    # STEP 2 -- Set up final DREAMSettings object
    ###############################################
    ds.eqsys.n_re.setEceff(RunawayElectrons.COLLQTY_ECEFF_MODE_FULL)
    if mode == MODE_FLUID:
        ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_ULTRA_RELATIVISTIC
        ds.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_FLUID_HESSLOW)
        ds.eqsys.n_re.setDreicer(RunawayElectrons.DREICER_RATE_NEURAL_NETWORK)

        # Include fluid hot-tail
        if not np.isscalar(T0r):
            rn0 = T0r
            n0 = nfree*np.ones(T0r.shape)
        else:
            rn0 = 0
            n0 = nfree

        nre = 0
        if nre0 is not None:
            nre = nre0

        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=n0-nre, rT0=T0r, T0=T0)
        ds.eqsys.n_re.setHottail(RunawayElectrons.HOTTAIL_MODE_ANALYTIC_ALT_PC)
        ds.eqsys.T_cold.setInitialProfile(T0, radius=T0r)
    else:
        if mode == MODE_ISOTROPIC:
            ds.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_FLUID_HESSLOW)
            ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
            ds.eqsys.T_cold.setInitialProfile(temperature=1)
        else:
            ds.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.01)
            ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
            ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_EXPLICIT, shape=FHot.PARTICLE_SOURCE_SHAPE_DELTA)

            ignorelist = ['n_i', 'N_i', 'W_i']

        ds.eqsys.n_re.setDreicer(RunawayElectrons.DREICER_RATE_DISABLED)
        
        # Kinetic ionization
        ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)
        #ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_FLUID)

    # Set self-consistent E and T evolution
    ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    ds.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT)
    if mode == MODE_KINETIC:
        ds.eqsys.E_field.setInitialProfile(Einit, radius=Einit_r)
    else:
        ds.eqsys.j_ohm.setInitialProfile(j0(j0r), radius=j0r, Ip0=Ip0)

    if tauwall is None:
        Vloop_wall_R0 = 0
        if Vloop is not None:
            Vloop_wall_R0 = Vloop / R0

        ds.eqsys.E_field.setBoundaryCondition(
            ElectricField.BC_TYPE_PRESCRIBED,
            V_loop_wall_R0=Vloop_wall_R0)
    elif Vloop is None or Vloop == 0:
        ds.eqsys.E_field.setBoundaryCondition(
            ElectricField.BC_TYPE_SELFCONSISTENT,
            inverse_wall_time=1/tauwall, R0=R0)
    else:
        ds.eqsys.E_field.setBoundaryCondition(
            ElectricField.BC_TYPE_TRANSFORMER,
            inverse_wall_time=1/tauwall, R0=R0,
            V_loop_wall_R0=Vloop/R0)

    # Runaway electron grid
    if withfre and not ds.runawaygrid.enabled:
        setup_runawaygrid(ds, equilibrium)

        if (nre0 is not None) and (nre0 != 0):
            #ds.eqsys.f_re.setInitType(FRe.INIT_AVALANCHE)
            #ds.eqsys.f_re.setInitType(FRe.INIT_FORWARD)
            ds.eqsys.f_re.setInitType(FRe.INIT_XI_NEGATIVE)

    MAXIMUM_IGNORABLE_ELECTRON_DENSITY = 1e5
    ds.solver.tolerance.set(reltol=reltol)
    ds.solver.tolerance.set(unknown='n_re', reltol=reltol, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
    ds.solver.tolerance.set(unknown='j_re', reltol=reltol, abstol=e*c*MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
    ds.solver.tolerance.set(unknown='f_re', reltol=reltol, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)

    ds.solver.setMaxIterations(500)

    if mode == MODE_ISOTROPIC:
        ds.solver.tolerance.set(unknown='f_hot', reltol=reltol, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
        ds.solver.tolerance.set(unknown='n_hot', reltol=reltol, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
        ds.solver.tolerance.set(unknown='j_hot', reltol=reltol, abstol=e*c*MAXIMUM_IGNORABLE_ELECTRON_DENSITY)

    ds.output.setTiming(True, True)
    
    if mode == MODE_KINETIC:
        ds.fromOutput(INITFILE, ignore=ignorelist)

    return ds


def simulate(ds1, mode, impurities, t_sim, dt0, dtmax, Drr=0,
    nre0=None, nre0_r=0,
    verboseIoniz=False, runIoniz=True, verboseTQ=False, runTQ=True,
    prefix='output/generic', extension='', **kwargs):
    """
    Run the disruption simulation (given a baseline simulation).

    :param impurities: List of impurity specifications (dict with keys: name, Z, n)
    """
    setname = lambda phase : longsetname(mode=mode, phase=phase, prefix=prefix, extension=extension)
    outname = lambda phase : longoutname(mode=mode, phase=phase, prefix=prefix, extension=extension)

    Ti0 = 1

    ##########################################
    # 1. Ionization phase of TQ (~1 µs)
    ##########################################
    # Add injected impurities
    for i in impurities:
        print('i[n]:', i['n'])
        ds1.eqsys.n_i.addIon(i['name'], Z=i['Z'], iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=i['n'], T=Ti0)

    # Prescribe heat diffusion?
    if Drr > 0:
        ds1.eqsys.T_cold.transport.prescribeDiffusion(drr=Drr)

    #ds1.timestep.setTmax(t_ioniz)
    #ds1.timestep.setNt(nt_ioniz)

    #if nt_ioniz > 5000:
    #    ds1.timestep.setNumberOfSaveSteps(5000)
    #ds1.timestep.setIonization(dtmax=5e-7, tmax=t_TQ, safetyfactor=800)
    ds1.timestep.setIonization(dt0=dt0, dtmax=dtmax, tmax=t_sim)

    ds1.solver.setVerbose(verboseIoniz)

    IONIZOUT = outname('1')
    ds1.save(setname('1'))

    if runIoniz:
        do1 = runiface(ds1, IONIZOUT)
    else:
        try: do1 = DREAMOutput(IONIZOUT)
        except: do1 = None


    ds2 = DREAMSettings(ds1)
    dBB = 4e-4 # Remnant heat transport after TQ
    ds2.eqsys.T_cold.transport.setMagneticPerturbation(dBB=dBB)

    ds2.timestep.setTmax(3e-2 - t_sim)
    ds2.timestep.setDt(dtmax)
    ds2.timestep.setType(TimeStepper.TYPE_CONSTANT)

    do2 = runiface(ds2, f'output_CQ.h5')

    """
    ##########################################
    # 2. Thermal and current quench
    ##########################################
    ds2 = DREAMSettings(ds1, keepignore=False)
    ds2.fromOutput(IONIZOUT)

    ds2.timestep.setTmax(t_TQ - t_ioniz)
    ds2.timestep.setNt(nt_TQ)

    if nt_TQ > 5000:
        ds2.timestep.setNumberOfSaveSteps(5000)
    else:
        ds2.timestep.setNumberOfSaveSteps(0)

    ds2.solver.setVerbose(verboseTQ)
    TQOUT = outname('2')
    ds2.output.setFilename(TQOUT)
    ds2.save(setname('2'))

    #ds2.solver.setDebug(savesystem=True, timestep=9, iteration=0)

    if runTQ:
        do2 = runiface(ds2, TQOUT)
    else:
        try: do2 = DREAMOutput(TQOUT)
        except: do2 = None

    return do1, do2, ds2
    """
    return do1, ds1


def plateau_simulation(dso, mode, tmax, index=0, **kwargs):
    """
    Continue the simulation through the plateau. This function can be
    called repeatedly.
    """
    setname = lambda phase : longsetname(mode=mode, phase=phase, prefix=kwargs['prefix'], extension=kwargs['extension'])
    outname = lambda phase : longoutname(mode=mode, phase=phase, prefix=kwargs['prefix'], extension=kwargs['extension'])

    ds = DREAMSettings(dso, keepignore=False)

    ds.timestep.setType(TimeStepper.TYPE_CONSTANT)
    ds.timestep.setTmax(tmax)
    ds.timestep.setNt(None)
    ds.timestep.setDt(1e-5)
    ds.timestep.setNumberOfSaveSteps(0)

    OUT = outname(f'{index+2:d}')
    ds.solver.setVerbose(kwargs['verbosePlateau'])
    ds.output.setFilename(OUT)
    ds.save(setname(f'{index+2:d}'))

    if kwargs['runPlateau']:
        do = runiface(ds, OUT)
    else:
        try: do = DREAMOutput(OUT)
        except: do = None

    return do


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Run a TCV disruption simulation")

    parser.add_argument('-c', '--cylindrical', help="Run in cylindrical geometry", dest="cylindrical", action="store_true")
    parser.add_argument('-e', '--extension', help="Append the given string to the end of all file names", dest="extension", action="store", type=str)
    parser.add_argument('-f', '--with-fre', help="Include a runaway electron distribution function", dest="withfre", action="store_true")
    parser.add_argument('--fluid', help="Run a simulation in pure fluid mode", dest="mode", action="store_const", const=MODE_FLUID)
    parser.add_argument('--isotropic', help="Run a simulation in isotropic mode", dest="mode", action="store_const", const=MODE_ISOTROPIC)
    parser.add_argument('--kinetic', help="Run a simulation in the fully kinetic mode", dest="mode", action="store_const", const=MODE_KINETIC)
    parser.add_argument('-r', '--run-from', help="Determines which simulation to start running from (0 = current, 1 = ionization, 2 = CQ)", dest="runfrom", action="store", type=int)
    parser.add_argument('-v', '--verbose', help="Select which simulations to make verbose (0 = current, 1 = ionization, 2 = CQ)", dest="verbose", nargs='*', type=int)

    parser.add_argument('--scan-current', help="Run a current density profile scan", dest="scan", action="store_const", const=SCAN_CURRENT)
    parser.add_argument('--scan-d2', help="Run a D2 injection scan", dest="scan", action="store_const", const=SCAN_D2)
    parser.add_argument('--scan-drr', help="Run a heat transport scan", dest="scan", action="store_const", const=SCAN_DRR)
    parser.add_argument('--scan-neon', help="Run a neon scan", dest="scan", action="store_const", const=SCAN_NEON)
    parser.add_argument('--scan-nre', help="Run a initial runaway density scan", dest="scan", action="store_const", const=SCAN_NRE)
    parser.add_argument('--scan-temperature', help="Run a central temperature scan", dest="scan", action="store_const", const=SCAN_TEMPERATURE)
    parser.add_argument('--scan-temperature-profile', help="Run a temperature profile scan", dest="scan", action="store_const", const=SCAN_TEMPERATURE_PROFILE)
    parser.add_argument('--scan-vloop', help="Run a scan in the edge loop voltage", dest="scan", action="store_const", const=SCAN_VLOOP)

    parser.set_defaults(cylindrical=False, mode=MODE_FLUID, scan=SCAN_NONE, runfrom=0, verbose=[], extension='', withfre=False)

    return parser.parse_args(args)

