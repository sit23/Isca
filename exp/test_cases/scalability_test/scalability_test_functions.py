import xarray as xar
import sys
import os
from isca import Experiment, IscaCodeBase, GreyCodeBase, FailedRunError, GFDL_BASE, DiagTable
from isca.util import exp_progress
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import json
import pdb

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/trip_test/'))
from trip_test_functions import get_nml_diag, define_simple_diag_table

def     run_scalability_test(test_case_name='frierson', list_of_num_procs=[4,8], resolutions_to_check=['T42', 'T85'], num_days=3, repo_to_use='git@github.com:execlim/Isca', commit_id ='HEAD'):

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
    total_time_short_arr = {}

    for resolution in resolutions_to_check:
        delta_t_arr_1={}
        delta_t_arr_2={}    
        compute_time_per_day_arr={}
        other_time = {}
        total_time = {}
        total_time_short = {}        

        for num_procs in list_of_num_procs:
            exp_name = test_case_name+'_scalability_test_'+str(num_procs)+'_'+resolution
            exp = Experiment(exp_name, codebase=cb)
            exp.namelist = nml_use.copy()
            exp.diag_table = diag_use
            exp.inputfiles = input_files_use

            exp.set_resolution(resolution)

            exp.update_namelist({
            'fms_nml': {
            'domains_stack_size':2000000,
            }})

            #Run first 'month' for one day, to make sure run that has its time measured has restarting involved.
            exp.update_namelist({
            'main_nml': {
            'days':1,
            }})


            try:
                # run with a progress bar
                with exp_progress(exp, description=str(num_procs)) as pbar:
                    exp.run(1, use_restart=False, num_cores=num_procs, overwrite_data=True)
                    
                exp.update_namelist({
                'main_nml': {
                'days':num_days,
                }})                    
                
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
            total_time_short[num_procs] = delta_t_arr_1[num_procs]/ exp.namelist['main_nml']['days']
                
        res_dict_per_day[resolution] = compute_time_per_day_arr
        restart_file_write_time[resolution] = other_time
        total_time_arr[resolution] = total_time
        total_time_short_arr[resolution] = total_time_short

    return res_dict_per_day, restart_file_write_time, total_time_arr, total_time_short_arr, commit_id

def create_output_folder(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used):

    GFDL_BASE_DIR = os.environ['GFDL_BASE']
    
    list_procs_str=[str(num) for num in list_of_num_procs]
    num_procs_str = '_'.join(list_procs_str)
    num_days=str(num_days)+'_days'
    res_str = '_'.join(resolutions_to_check)
    
    
    extra_dir = 'output_plots/'+test_case_name+'/'+'procs_'+num_procs_str+'_'+res_str+'_'+num_days+'_commit_'+commit_id_used
    directory = GFDL_BASE_DIR + '/exp/test_cases/scalability_test/'+extra_dir

    if not os.path.exists(directory):
        os.makedirs(directory)

    return directory

def output_to_file(out_dir,linr_compute_time, linr_write_time, total_time, total_time_half_length):
    
    with open(out_dir+'/values_linr_compute_time.json', 'w') as f:  json.dump(linr_compute_time, f)
    with open(out_dir+'/values_linr_write_time.json', 'w') as f: json.dump(linr_write_time, f)    
    with open(out_dir+'/values_total_time.json', 'w') as f: json.dump(total_time, f)    
    with open(out_dir+'/values_total_time_half_length.json', 'w') as f: json.dump(total_time_half_length, f) 

def attempt_read_in(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used):  

    out_dir = create_output_folder(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used)

    with open(out_dir+'/values_linr_compute_time.json', 'r') as f: linr_compute_time=json.load(f)
    with open(out_dir+'/values_linr_write_time.json', 'r') as f: linr_write_time=json.load(f)
    with open(out_dir+'/values_total_time.json', 'r') as f: total_time = json.load(f)
    with open(out_dir+'/values_total_time_half_length.json', 'r') as f: total_time_half_length = json.load(f)
    
    make_keys_numbers(linr_compute_time)
    make_keys_numbers(linr_write_time)
    make_keys_numbers(total_time)
    make_keys_numbers(total_time_half_length)
    
    return linr_compute_time, linr_write_time, total_time, total_time_half_length, commit_id_used

def make_keys_numbers(input_dict):
    for res in input_dict.keys():
        input_dict[res] = {int(key):input_dict[res][key] for key in input_dict[res].keys()}

def create_output(linr_compute_time, linr_write_time, total_time, total_time_half_length, test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used):
    
    out_dir = create_output_folder(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used)

    output_to_file(out_dir,linr_compute_time, linr_write_time, total_time, total_time_half_length)
    
    print('Creating plots')
    
    for resolution in resolutions_to_check:
        time_dict = linr_compute_time[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution, marker='x')

    theoretical_perfect = [entry/x_values[0] for entry in x_values]
    plt.plot(x_values, theoretical_perfect, label='Perfect scaling', marker='x')
    plt.title('Time per iteration for '+test_case_name+' including output writing\nbut excluding restart-file time (linr calc)')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with\nfewest number of cores')

    plt.legend()

    plt.savefig(out_dir+'/linr_times_per_iteration.pdf')

    plt.figure()
    for resolution in resolutions_to_check:
        time_dict = total_time[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution, marker='x')

    theoretical_perfect = [entry/x_values[0] for entry in x_values]
    plt.plot(x_values, theoretical_perfect, label='Perfect scaling', marker='x')
    plt.title('Time per iteration for '+test_case_name+'\nincluding output writing and restart-file time')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with\nfewest number of cores')

    plt.legend()

    plt.savefig(out_dir+'/total_times_per_iteration.pdf')

    plt.figure()
    for resolution in resolutions_to_check:
        time_dict = total_time_half_length[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution, marker='x')

    theoretical_perfect = [entry/x_values[0] for entry in x_values]
    plt.plot(x_values, theoretical_perfect, label='Perfect scaling', marker='x')
    plt.title('Time per iteration for '+test_case_name+'\nincluding output writing and restart-file time')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with\nfewest number of cores')

    plt.legend()

    plt.savefig(out_dir+'/total_times_short_per_iteration.pdf')

    plt.figure()
    for resolution in resolutions_to_check:
        time_dict = linr_write_time[resolution]
        x_values = list(time_dict.keys())
        y_values = list(time_dict.values())
    
        y_values_scaled = [y_values[0]/entry for entry in y_values]
    
        plt.plot(x_values, y_values_scaled, label=resolution, marker='x')
    plt.legend()
    plt.title('Day-number independent part of run time for '+test_case_name+' (linr calc)')
    plt.xlabel('# of cores')
    plt.ylabel('Number of times faster than case with\nfewest number of cores')

    plt.savefig(out_dir+'/linr_times_that_dont_scale_with_day_number.pdf')

def scalability_control(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id='HEAD', force_rerun=False):

    if force_rerun:
        linr_compute_time, linr_write_time, total_time, total_time_half_length, commit_id_used= run_scalability_test(test_case_name, list_of_num_procs, resolutions_to_check, num_days, commit_id=commit_id)    
    else:    
        try:
            print('Attempting to read previous values')
            linr_compute_time, linr_write_time, total_time, total_time_half_length, commit_id_used = attempt_read_in(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id)
            print('Successfully read previous values')        

        except:
            print('Cannot read previous values - running tests')    
            linr_compute_time, linr_write_time, total_time, total_time_half_length, commit_id_used= run_scalability_test(test_case_name, list_of_num_procs, resolutions_to_check, num_days, commit_id=commit_id)

    create_output(linr_compute_time, linr_write_time, total_time, total_time_half_length, test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id_used)


if __name__=="__main__":
    test_case_name_list=['held_suarez', 'realistic_continents_variable_qflux']
    list_of_num_procs = [4, 8, 16,32]
    num_days=15
    resolutions_to_check = ['T42', 'T85', 'T170']
    commit_id = '5868cfd'
    force_rerun = False

    for test_case_name in test_case_name_list:
        scalability_control(test_case_name, list_of_num_procs, num_days, resolutions_to_check, commit_id=commit_id, force_rerun=force_rerun)