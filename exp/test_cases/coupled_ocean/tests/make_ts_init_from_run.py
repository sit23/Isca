"""Initial ocean temperature/salinity for INIT_LAYERS_FROM_Z_FILE from a MOM6 run's state.

Takes the zonal-mean potential temperature of EXP's month N (ocean_month_z.nc, on z
levels) and the mean sea-surface salinity, puts them on a global 1-degree grid (as
mom6.write_initial_ts does) with latitudes beyond the run's domain given the values of
its last row, and writes OUT. Used to start runs on other ocean grids (e.g. tripolar)
from the same ocean state.

Usage: python make_ts_init_from_run.py EXP N OUT
"""
import os
import sys

import netCDF4
import numpy as np
import xarray as xr

D = os.environ['GFDL_DATA']
exp, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
odir = os.path.join(D, exp, 'run%04d' % n, 'ocean')
z = xr.open_dataset(os.path.join(odir, 'ocean_month_z.nc'), decode_times=False).isel(time=0)
m = xr.open_dataset(os.path.join(odir, 'ocean_month.nc'), decode_times=False).isel(time=0)
st = xr.open_dataset(os.path.join(odir, 'ocean_static.nc'))

lat_run = st.geolat.mean('xh').values
temp_run = z.thetao.mean('xh').values            # (z, lat)
depth = z.z_l.values
salt = float(np.nansum(m.sos * st.area_t) / np.nansum(st.area_t.where(np.isfinite(m.sos))))

lat = np.arange(-89.5, 90., 1.)
lon = np.arange(0.5, 360., 1.)
temp = np.empty((len(depth), len(lat)))
for k in range(len(depth)):
    col = temp_run[k]
    ok = np.isfinite(col)
    temp[k] = np.interp(lat, lat_run[ok], col[ok])   # constant beyond the run's last rows
# below the deepest valid level, keep the level above
for k in range(1, len(depth)):
    bad = ~np.isfinite(temp[k])
    temp[k, bad] = temp[k - 1, bad]

with netCDF4.Dataset(out, 'w') as nc:
    nc.createDimension('lon', len(lon)); nc.createDimension('lat', len(lat))
    nc.createDimension('depth', len(depth))
    for name, data, units, axis in (('lon', lon, 'degrees_east', 'X'), ('lat', lat, 'degrees_north', 'Y'),
                                    ('depth', depth, 'm', 'Z')):
        v = nc.createVariable(name, 'f8', (name,))
        v.units = units; v.axis = axis
        if axis == 'Z':
            v.positive = 'down'
        v[:] = data
    t = nc.createVariable('ptemp', 'f8', ('depth', 'lat', 'lon'), fill_value=-1.e20)
    t.missing_value = -1.e20; t.units = 'degC'
    t[:] = np.broadcast_to(temp[:, :, None], (len(depth), len(lat), len(lon)))
    s = nc.createVariable('salt', 'f8', ('depth', 'lat', 'lon'), fill_value=-1.e20)
    s.missing_value = -1.e20; s.units = 'psu'
    s[:] = salt
print('wrote', out, 'salinity %.3f, surface T range %.2f..%.2f' % (salt, np.nanmin(temp[0]), np.nanmax(temp[0])))
