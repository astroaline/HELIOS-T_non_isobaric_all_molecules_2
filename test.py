import numpy as np
from input import *
from asthetics import *
import data_setup
import ns_setup
import model
from matplotlib import pyplot as plt
import cornerplot
import time
import json
import os
import csv
import pdb

start = time.time()

n_params = len(parameters)


x, x_full, opacity_grid, bin_indices, ydata, yerr, wavelength_centre, wavelength_err = data_setup.data()
len_x = len(x)

new_param_dict = parameter_dict
new_param_dict['T'] = 1000
new_param_dict['log_xh2o'] = -2
new_param_dict['log_kappa_cloud'] = -5


## compute model for retrieved results ##
y = model.Model(len_x, x_full, bin_indices, new_param_dict, opacity_grid)
yfit = y.transit_depth()
yfit_binned = y.binned_model()


