import numpy as np
import matplotlib.pyplot as plt
import pdb

def convert_zhalf(zhalf, npfull, p_surf, add_extras=True,):

    phalf_array = p_surf*np.exp(-zhalf)
    if add_extras:
        phalf_array = np.concatenate((phalf_array, np.array([0.])))            
        zhalf_array = np.concatenate((zhalf, np.array([np.inf])))  
    else:
        zhalf_array=zhalf

    pfull_array = np.asarray([0.5*(phalf_array[pidx]+phalf_array[pidx-1]) for pidx in range(1,npfull+1)])

    zfull_array = np.asarray([0.5*(zhalf_array[pidx]+zhalf_array[pidx-1]) for pidx in range(1,npfull+1)])

    return phalf_array, pfull_array, zhalf_array, zfull_array

def calc_levels(npfull, p_surf):

    zhalf_array = np.linspace(0., 5., num=npfull, endpoint=True) #only use npfull so that 0. can be added at the end to make nphalf
   
    phalf_array, pfull_array, zhalf_array, zfull_array = convert_zhalf(zhalf_array, npfull, p_surf)

    return pfull_array, phalf_array, zfull_array, zhalf_array

def calc_merged_levels(n_low, n_high, z_min=1., z_max=2., p_surf=12000.):

    pfull_array, phalf_array, zfull_array, zhalf_array = calc_levels(n_low, p_surf)

    pfull_array_high, phalf_array_high, zfull_array_high, zhalf_array_high = calc_levels(n_high, p_surf)

    reduced_zhalf_array = zhalf_array[np.where(np.logical_or(zhalf_array>z_max, zhalf_array<z_min,))]

    required_zhalf_array_high = zhalf_array_high[np.where(np.logical_and(zhalf_array_high<=z_max, zhalf_array_high>=z_min,))]

    merged_zhalf_array = np.concatenate([reduced_zhalf_array, required_zhalf_array_high])

    sorted_zhalf_array = np.sort(merged_zhalf_array)

    plt.figure()
    plt.plot(zhalf_array, zhalf_array, marker='+', linestyle='none')
    plt.plot(zhalf_array_high, zhalf_array_high, marker='x', linestyle='none')
    plt.plot(sorted_zhalf_array, sorted_zhalf_array, marker='o', linestyle='none')

    npfull_sorted = (sorted_zhalf_array.shape[0]-1)

    phalf_array_conv, pfull_array_conv, zhalf_array_conv, zfull_array_conv = convert_zhalf(sorted_zhalf_array, npfull_sorted, p_surf, add_extras=False)

    return phalf_array_conv, pfull_array_conv, zhalf_array_conv, zfull_array_conv


if __name__=="__main__":
    phalf_array_conv, pfull_array_conv, zhalf_array_conv, zfull_array_conv = calc_merged_levels(21, 61, p_surf=10000.)    