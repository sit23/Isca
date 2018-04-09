import xarray as xar
import sys
import os
from isca import Experiment, IscaCodeBase, GreyCodeBase, FailedRunError, GFDL_BASE, DiagTable
from isca.util import exp_progress
import matplotlib.pyplot as plt
import time
import pdb

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/trip_test/'))
from trip_test_functions import get_nml_diag, define_simple_diag_table

test_case_name='frierson'
list_of_num_procs = [4, 8,16,32]
num_months=2
repo_to_use='git@github.com:execlim/Isca'
resolutions_to_check = ['T42', 'T85']
do_output = True

nml_use, input_files_use  = get_nml_diag(test_case_name)
if do_output:
    diag_use = define_simple_diag_table()
    exp_name_diag=''
else:
    diag_use = DiagTable()
    exp_name_diag='_no_output'
    
test_pass = True
run_complete = True

#Do the run for each of the commits in turn
s='5868cfd'
cb = GreyCodeBase(repo=repo_to_use, commit=s)
cb.compile()

res_dict={}

for resolution in resolutions_to_check:
    delta_t_arr={}

    for num_procs in list_of_num_procs:
        exp_name = test_case_name+'_scalability_test_'+str(num_procs)+'_'+resolution+exp_name_diag
        exp = Experiment(exp_name, codebase=cb)
        exp.namelist = nml_use.copy()
        exp.diag_table = diag_use
        exp.inputfiles = input_files_use

        exp.set_resolution(resolution)

        #Only run for 3 days to keep things short.
        exp.update_namelist({
        'main_nml': {
        'days': 3,
        }})

        try:
            # run with a progress bar
            with exp_progress(exp, description=str(num_procs)) as pbar:
                exp.run(1, use_restart=False, num_cores=num_procs, overwrite_data=True)
            start_time = time.time()
            for i in range(2,num_months+1):
                with exp_progress(exp, description=str(num_procs)) as pbar:        
                    ste = exp.run(i, num_cores=num_procs, overwrite_data=True)                                    
        except FailedRunError as e:
            #If run fails then test automatically fails
            run_complete = False
            test_pass = False
            continue
        end_time = time.time()
        delta_t_arr[num_procs] = (end_time-start_time)/(((num_months+1)-2)* exp.namelist['main_nml']['days'])
    res_dict[resolution] = delta_t_arr
    
    
for resolution in resolutions_to_check:
    time_dict = res_dict[resolution]
    x_values = list(time_dict.keys())
    y_values = list(time_dict.values())
    
    y_values_scaled = [y_values[0]/entry for entry in y_values]
    
    plt.plot(x_values, y_values_scaled, label=resolution)

plt.legend()