import numpy as np
import os

## Constants ##

kboltz = 1.38064852e-16    # Boltzmann's constant
amu = 1.660539040e-24      # atomic mass unit
gamma = 0.57721
rjup = 7.1492e9            # equatorial radius of Jupiter
rsun = 6.9566e10           # solar radius
rearth = 6.378e8            # earth radius
pressure_probed = 1e-2      # probed pressure in bars
# pressure_cia = 1e-2         # pressure for cia in bars
# m = 2.4*amu                 # assummed hydrogen-dominated atmosphere
m_water = 18.0*amu          # mean molecular mass of any molecules you want to consider
m_cyanide = 27.0*amu
m_ammonia = 17.0*amu
m_methane = 16.0*amu
m_carbon_monoxide = 28.0*amu


## Planet Data ##

planet_name = 'WASP-12b'

g = 977
g_uncertainty = 67
rstar = 1.57
rstar_uncertainty = 0.07
r0 = 1.75    # Rp = [1.70,1.88]
r0_uncertainty = 0.15    # R0 = [1.60,1.90]

wavelength_bins = np.array([0.838, 0.896, 0.954, 1.012, 1.070, 1.112, 1.182, 1.251, 1.320, 1.389, 1.458, 1.527, 1.597, 1.666])
transit_depth = np.array([1.4441, 1.4422, 1.4402, 1.4428, 1.4391, 1.4386, 1.4365, 1.4327, 1.4582, 1.4600, 1.4530, 1.4475, 1.4332])
transit_depth_error = np.array([0.0069, 0.0055, 0.0052, 0.0051, 0.0053, 0.0047, 0.0045, 0.0041, 0.0040, 0.0043, 0.0045, 0.0058, 0.0055])

pmin = 1e-6



## Retrieval info ##

model_name = 'greycloud'

molecules = ['1H2-16O__POKAZATEL_e2b', '12C-1H4__YT10to10_e2b']  # list of molecules (determines which opacity tables are loaded)
parameters = ["T", "log_xh2o", "log_xch4", "log_kappa_cloud", "R0", "Rstar", "G"]   # parameters you wish to retrieve (MUST MATCH MOLECULES)
res = 2         # resolution used for opacities
live = 1000     # live points used in nested sampling
wavenumber=True     # True if opacity given in terms of wavenumber, False if wavelength

priors = {"T": [2700, 200], "log_xh2o": [13,-13], "log_xch4": [13,-13], "log_xco": [13,-13], "log_kappa_cloud": [14,-12],
          "log_P0": [4,-1], "log_kappa_0": [9,-10], "Q0": [99,1], "a": [3,3], "R0": [2*r0_uncertainty, r0-r0_uncertainty],
          "log_r_c": [6,-7], "log_p_cia": [3,-3], "Rstar": [2*rstar_uncertainty,rstar-rstar_uncertainty], "G": [2*g_uncertainty,g-g_uncertainty], "line": [5,0]} # priors for all possible parameters



## info for all possible parameters ##
molecular_name_dict = {'1H2-16O__POKAZATEL_e2b': 'water', '12C-1H4__YT10to10_e2b': 'methane', '12C-16O__HITEMP2010_e2b': 'carbon_monoxide'}  # dictionary list of all possible molecules and corresponding names
molecular_abundance_dict = {'1H2-16O__POKAZATEL_e2b': 'log_xh2o', '12C-1H4__YT10to10_e2b': 'log_xch4', '12C-16O__HITEMP2010_e2b': 'log_xco'}  # dictionary list of all possible molecules and corresponding abundance names

parameter_dict = {"T": 1000, "log_xh2o": "Off", "log_xch4": "Off", "log_xco": "Off", "log_kappa_cloud": "Off", "R0": r0, "Rstar": rstar,
                  "log_P0": 1, "log_kappa_0": "Off", "Q0": "Off", "a": "Off", "log_r_c": "Off", "log_p_cia": -2, "G": g, "line": "Off"}    # default parameter values used if not retrieved

molecular_mass_dict = {'1H2-16O__POKAZATEL_e2b': m_water, '12C-1H4__YT10to10_e2b': m_methane, '12C-16O__HITEMP2010_e2b': m_carbon_monoxide}   # dictionary of molecules and their mean molecular masses
temperature_array = np.r_[50:700:50, 700:1500:100, 1500:3100:200]
# temperature_array = np.array([1500, 1700, 1900, 2100])
# temp_dict = {'01': temperature_array[9:], '12C-1H4__YT10to10_e2b': temperature_array[9:], '12C-16O__HITEMP2010_e2b': temperature_array}   # temperature values for corresponding opacity tables
temperature_array_cia = np.r_[200:3025:25]          # temperature array for CIA table
opacity_path = os.environ['HOME'] + "/Desktop/PhD/OPACITIES/"  # path to opacity binary files
cia_path = os.environ['HOME'] + "/Desktop/PhD/HITRAN/"      # path to CIA files
