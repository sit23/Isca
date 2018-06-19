import numpy as np
import matplotlib.pyplot as plt

ice_albedo_value = 0.7
albedo_value = 0.31

for ice_albedo_feedback_width in [1., 2., 5.,10.,20.]:

    t_surf = np.arange(100.,350,1.)

    ice_concentration = 0.5*(np.tanh(-1.*(t_surf-273.)/ice_albedo_feedback_width)+1.0)

    albedo = albedo_value + ice_concentration*(ice_albedo_value - albedo_value)

    plt.figure(1)
    plt.plot(t_surf, ice_concentration, label = ice_albedo_feedback_width)
    
    plt.figure(2)
    plt.plot(t_surf, albedo, label = ice_albedo_feedback_width)
    
plt.figure(1)    
plt.legend()
plt.figure(2)
plt.legend()

plt.show()
