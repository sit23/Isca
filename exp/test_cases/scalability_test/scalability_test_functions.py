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

def     run_scalability_test(test_case_name='frierson', list_of_num_procs=[4,8], resolutions_to_check=['T42', 'T85'], num_days=3, repo_to_use='git@github.com:execlim/Isca', commit_id ='5868cfd'):

    nml_use, input_files_use  = get_nml_diag(test_case_name)
    diag_use = define_simple_diag_table()
    
    test_pass = True
    run_complete = True

    #Do the run for each of the commits in turn
    cb = GreyCodeBase(repo=repo_to_use, commit=commit_id)
    cb.compile()

    res_dict_per_day={}
    restart_file_write_time = {}
    total_time_arr = {}

    for resolution in resolutions_to_check:
        delta_t_arr_1={}
        delta_t_arr_2={}    
        compute_time_per_day_arr={}
        other_time = {}
        total_time = {}
    


        for num_procs in list_of_num_procs:
            exp_name = test_case_name+'_scalability_test_'+str(num_procs)+'_'+resolution
            exp = Experiment(exp_name, codebase=cb)
            exp.namelist = nml_use.copy()
            exp.diag_table = diag_use
            exp.inputfiles = input_files_use

            exp.set_resolution(resolution)

            #Only run for 3 days to keep things short.
            exp.update_namelist({
            'main_nml': {
            'days':num_days,
            }})

            try:
                # run with a progress bar
                with exp_progress(exp, description=str(num_procs)) as pbar:
                    exp.run(1, use_restart=False, num_cores=num_procs, overwrite_data=True)
                start_time = time.time()
                with exp_progress(exp, description=str(num_procs)) as pbar:        
                    ste = exp.run(2, num_cores=num_procs, overwrite_data=True)                                    
                end_time_1_iteration = time.time()
                exp.update_namelist({
                'main_nml': {
                'days': num_days*2,
                }})
                start_time_2_iteration = time.time()
                with exp_progress(exp, description=str(num_procs)) as pbar:        
                    ste = exp.run(3, num_cores=num_procs, overwrite_data=True)                                                
                end_time_2_iteration = time.time()
            
            except FailedRunError as e:
                #If run fails then test automatically fails
                run_complete = False
                test_pass = False
                continue
            delta_t_arr_1[num_procs] = (end_time_1_iteration-start_time)
            delta_t_arr_2[num_procs] = (end_time_2_iteration-start_time_2_iteration)
        
            compute_time_per_day_arr[num_procs] = (delta_t_arr_2[num_procs]-delta_t_arr_1[num_procs]) / (0.5 * exp.namelist['main_nml']['days'])
        
            other_time[num_procs] = delta_t_arr_1[num_procs] - (compute_time_per_day_arr[num_procs] * (0.5 * exp.namelist['main_nml']['days']))
        
            total_time[num_procs] = delta_t_arr_2[num_procs]/ exp.namelist['main_nml']['days']
                
        res_dict_per_day[resolution] = compute_time_per_day_arr
        restart_file_write_time[resolution] = other_time
        total_time_arr[resolution] = total_time

    return res_dict_per_day, restart_file_write_time, total_time_arr

def create_output(linr_compute_time, linr_write_time, total_time):
    
    for resolution in resolutions_to_check:
        time_dict = res_dict_per_day[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution)

    theoretical_perfect = [entry/x_values[0] for entry in x_values]
    plt.plot(x_values, theoretical_perfect, label='Perfect scaling')
    plt.title('Time per iteration, including output writing but excluding restart-file time (linr calc)')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with fewest number of cores')

    plt.legend()

    plt.figure()
    for resolution in resolutions_to_check:
        time_dict = total_time_arr[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution)

    theoretical_perfect = [entry/x_values[0] for entry in x_values]
    plt.plot(x_values, theoretical_perfect, label='Perfect scaling')
    plt.title('Total Time per iteration, including output writing and restart-file time')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with fewest number of cores')

    plt.legend()

    plt.figure()
    for resolution in resolutions_to_check:
        time_dict = restart_file_write_time[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution)
    plt.legend()
    plt.title('Day-number independent part of run time (linr calc)')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with fewest number of cores')


if __name__=="__main__":
    test_case_name='held_suarez'
    list_of_num_procs = [4, 8, 16,32]
    num_days=3
    resolutions_to_check = ['T42', 'T85']

    linr_compute_time, linr_write_time, total_time = run_scalability_test(test_case_name, list_of_num_procs, resolutions_to_check, num_days)
    create_output(linr_compute_time, linr_write_time, total_time)
