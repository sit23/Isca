import numpy

def dry_adiabatic_lapse_rate(grav=9.81, specific_gas_constant_dry=287.04, kappa=2./7.):

    cp_air = specific_gas_constant_dry/kappa
    lapse_rate = grav/cp_air
    
    return lapse_rate*1000. # Multiply by 1000. so units are k/km.

def moist_adiabatic_lapse_rate(grav = 9.8, mixing_ratio = 1.e-3, specific_gas_constant_dry=287.04, kappa=2./7., temperature=283., specific_gas_constant_moist = 461.5):

    latent_heat_vaporization = 2.501e6
    cp_air = specific_gas_constant_dry/kappa

    lapse_rate_numerator = grav * (1. + (latent_heat_vaporization*mixing_ratio)/(specific_gas_constant_dry*temperature))
    lapse_rate_denominator = (cp_air + (latent_heat_vaporization**2.*mixing_ratio)/(specific_gas_constant_moist*temperature**2.))

    lapse_rate = lapse_rate_numerator / lapse_rate_denominator
    
    return lapse_rate*1000. # Multiply by 1000. so units are k/km.
    
if __name__=="__main__":

    #Earth
    dalr = dry_adiabatic_lapse_rate()
    malr = moist_adiabatic_lapse_rate()
    
    print('Planet    ', 'dry alr        ', 'moist alr       ')
    print('Earth     ', dalr, malr)
    
    #Jupiter
    dalr = dry_adiabatic_lapse_rate(grav=26.0, specific_gas_constant_dry=3605.38, kappa=2./7.)
    solar_mixing_ratio_wv = 1.7e-3
    for multiplier in [1.,2.,4.]:
        malr = moist_adiabatic_lapse_rate(grav = 26.0, mixing_ratio = solar_mixing_ratio_wv*multiplier, specific_gas_constant_dry=3605.38, kappa=2./7., temperature=220.)
        
        print('Jupiter_'+str(int(multiplier))+'s', dalr, malr)
    

