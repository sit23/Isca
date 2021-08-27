from netCDF4 import Dataset  
from plevel_fn import plevel_call, daily_average, join_files, two_daily_average, monthly_average
import sys
import os
import time
import pdb
import subprocess
import numpy as np

start_time=time.time()
#base_dir=os.environ['GFDL_DATA']
#base_dir = '/disca/share/sit204/data_from_isca_cpu/cssp_laura_exps/'
base_dir = '/disca/share/sit204/data_from_isca_cpu/cssp_perturb_exps/control/soc_ga3_do_simple_false_cmip_o3_bucket/'
#base_dir='/scratch/sit204/mounts/gv4/sit204/data_isca/'
#exp_name_list = ['soc_ga3_files_smooth_topo_fftw_mk1_fresh_compile_long', 'soc_ga3_files_smooth_topo_old_fft_mk2_long']
#exp_name_list = [f'soc_ga3_do_simple_false_cmip_o3_bucket_perturbed_ens_{f}' for f in range(200,400)]
#exp_name_list = ['soc_ga3_do_simple_false_cmip_o3_bucket_qflux_co2_400_mid_alb_eur_anom_mk2', 'soc_ga3_do_simple_false_cmip_o3_bucket_qflux_co2_400_mid_alb_asia_anom_mk2', 'soc_ga3_do_simple_false_cmip_o3_bucket_qflux_co2_400_mid_alb']
#exp_name_list = ['soc_ga3_do_simple_false_cmip_o3_bucket_qflux_co2_400_mid_alb_eur_anom_shifted_mk2']
#exp_name_list = ['control_exp']
##exp_name_list = ['soc_ga3_files_smooth_topo_mk1_fresh_compile_long_cmip5_ozone']
#exp_name_list = ['soc_ga3_files_smooth_topo_mk1_fresh_compile_long_ice_albedo_not_land_mask_ie_sst_jra_55_8_sbm_do_simple_false_cmip5_o3_bucket']
#exp_name_list = ['soc_ga3_files_smooth_topo_mk1_fresh_compile_long_ice_albedo_not_land_mask_ice_sst_jra_55_8_sbm_do_simple_false_cmip5_o3_bucket']
exp_name_list = ['control_exp']
avg_or_daily_list=['monthly']
start_file=3001
end_file=5500
nfiles=(end_file-start_file)+1

do_extra_averaging=False #If true, then 6hourly data is averaged into daily data using cdo
group_months_into_one_file=True # If true then monthly data files and daily data files are merged into one big netcdf file each.
overwrite_previous_combined_files=True
n_splits_of_combined_nc_file = 20
level_set='standard' #Default is the standard levels used previously. ssw_diagnostics are the ones blanca requested for MiMa validation
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 
all_vars=True
#var_names_list = 'slp height precipitation vcomp ucomp temp_2m temp div flux_t flux_lhe bucket_depth'
var_names_list = 'slp height temp_2m'

try:
    out_dir
except:
    out_dir = base_dir

plevs={}
var_names={}

if level_set=='standard':

    plevs['monthly']=' -p "3 16 51 138 324 676 1000 1266 2162 3407 5014 6957 9185 10000 11627 14210 16864 19534 20000 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 70000 73915 78095 82510 85000 87175 92104 97312"'

    plevs['timestep']=' -p "3 16 51 138 324 676 1000 1266 2162 3407 5014 6957 9185 10000 11627 14210 16864 19534 20000 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 70000 73915 78095 82510 85000 87175 92104 97312"'

    plevs['pentad']=' -p "3 16 51 138 324 676 1000 1266 2162 3407 5014 6957 9185 10000 11627 14210 16864 19534 20000 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 70000 73915 78095 82510 85000 87175 92104 97312"'

    plevs['6hourly']=' -p "1000 10000 25000 50000 85000 92500"'
    plevs['daily']  =' -p "1000 10000 25000 50000 85000 92500"'
    
    if all_vars:
        var_names['monthly']='-a slp height'
        var_names['pentad']='-a slp height'
    else:
        var_names['monthly']=var_names_list
        var_names['pentad']=var_names_list
    var_names['timestep']='-a'
    var_names['6hourly']='ucomp slp height vor t_surf vcomp omega'
    var_names['daily']='ucomp slp height vor t_surf vcomp omega temp'

    if all_vars:
        file_suffix='_interp_new_height_temp_not_below_ground'
    else:
        file_suffix='_interp_new_height_temp_not_below_ground_subset_vars_5'

elif level_set=='ssw_diagnostics':
    plevs['6hourly']=' -p "1000 10000"'
    var_names['monthly']='ucomp temp height'
    var_names['6hourly']='ucomp vcomp temp'
    file_suffix='_bl'
    
elif level_set=='tom_diagnostics':
    var_names['daily']='height temp'
    plevs['daily']=' -p "10 30 100 300 500 700 1000 3000 5000 7000 10000 15000 20000 25000 30000 40000 50000 60000 70000 75000 80000 85000 90000 95000 100000"'
    mask_below_surface_set='-x '
    file_suffix='_tom_mk2'


for exp_name in exp_name_list:
    for n in range(nfiles):
        for avg_or_daily in avg_or_daily_list:
            print(n+start_file)
            
            nc_file_in = base_dir+'/'+exp_name+'/run%04d'%(n+start_file)+'/atmos_'+avg_or_daily+'.nc'
            nc_file_out = out_dir+'/'+exp_name+'/run%04d'%(n+start_file)+'/atmos_'+avg_or_daily+file_suffix+'.nc'

            if not os.path.isfile(nc_file_out):
                plevel_call(nc_file_in,nc_file_out, var_names = var_names[avg_or_daily], p_levels = plevs[avg_or_daily], mask_below_surface_option=mask_below_surface_set)
            if do_extra_averaging and avg_or_daily=='6hourly':
                nc_file_out_daily = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_daily'+file_suffix+'.nc'
                daily_average(nc_file_out, nc_file_out_daily)
            if do_extra_averaging and avg_or_daily=='pentad':
                nc_file_out_daily = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_monthly'+file_suffix+'.nc'
                monthly_average(nc_file_out, nc_file_out_daily, adjust_time = True)                
#            if do_extra_averaging and avg_or_daily=='6hourly':
#                nc_file_out_two_daily = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_two_daily'+file_suffix+'.nc'
#                two_daily_average(nc_file_out, nc_file_out_two_daily, avg_or_daily)

if group_months_into_one_file:
    avg_or_daily_list_together=avg_or_daily_list


    for exp_name in exp_name_list:
        for avg_or_daily in avg_or_daily_list_together:
            
            for itr_idx in range(n_splits_of_combined_nc_file):
                nc_file_string=''
                n_files_in_each = int(np.ceil(nfiles / n_splits_of_combined_nc_file))

                start_file_itr = start_file + (itr_idx)*n_files_in_each
                end_file_itr = int(np.min([start_file + ((itr_idx+1)*n_files_in_each - 1), end_file]))

                for n in range(start_file_itr, end_file_itr+1):
                    nc_file_in = base_dir+'/'+exp_name+'/run%04d'%(n)+'/atmos_'+avg_or_daily+file_suffix+'.nc'
                    nc_file_string=nc_file_string+' '+nc_file_in
                nc_file_out=base_dir+'/'+exp_name+'/atmos_'+avg_or_daily+f'_together_{start_file_itr}_{end_file_itr}'+file_suffix+'.nc'
                if not os.path.isfile(nc_file_out) or overwrite_previous_combined_files:
                    join_files(nc_file_string,nc_file_out)

print('execution time', time.time()-start_time)



