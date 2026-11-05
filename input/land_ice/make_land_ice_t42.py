"""Ice-sheet fraction of the land on Isca's T42 grid, from ERA5 snow depth.

ERA5 holds snow depth on glaciers and ice sheets at 10 m water equivalent, so a 1-degree cell
whose snow depth stays above SD_MIN every month of 1991-2020 is ice sheet (Antarctica,
Greenland, the larger Arctic ice caps). Land at 1 degree is where ERA5 has no SST. For each
T42 cell: land_ice = (ice-sheet 1-degree cells) / (land 1-degree cells), by nearest T42 cell.

Usage: python make_land_ice_t42.py ERA5_SD.nc ERA5_SST.nc LAND_T42.nc OUT.nc
"""
import sys
import numpy as np
import xarray as xr

SD_MIN = 1.0   # m water equivalent
sd_f, sst_f, land_f, out_f = sys.argv[1:5]
sd = xr.open_dataset(sd_f).sd
sdmin = sd.groupby('valid_time.month').mean().min('month').values
sst = xr.open_dataset(sst_f).sst.isel(valid_time=0).values
lat1, lon1 = sd.latitude.values, sd.longitude.values
land1 = ~np.isfinite(sst)
ice1 = land1 & (sdmin >= SD_MIN)
w1 = np.cos(np.radians(lat1))[:, None] * np.ones(len(lon1))

t42 = xr.open_dataset(land_f)
lat, lon = t42.lat.values, t42.lon.values
jj = np.abs(lat1[:, None] - lat[None, :]).argmin(1)
dl = np.abs(((lon1[:, None] - lon[None, :]) + 180.) % 360. - 180.)
ii = dl.argmin(1)
J, I = np.meshgrid(jj, ii, indexing='ij')
nland = np.zeros((len(lat), len(lon)))
nice = np.zeros_like(nland)
np.add.at(nland, (J, I), w1 * land1)
np.add.at(nice, (J, I), w1 * ice1)
frac = np.where(nland > 0, nice / np.maximum(nland, 1e-12), 0.)
frac = np.where(t42.land_mask.values > 0, frac, 0.)
ds = xr.Dataset({'land_ice': (('lat', 'lon'), frac.astype('f4'))}, coords={'lat': lat, 'lon': lon})
ds.land_ice.attrs = dict(long_name='fraction of the land covered by ice sheet', units='1',
                         source='ERA5 snow depth >= %g m w.e. every month, 1991-2020' % SD_MIN)
ds.to_netcdf(out_f)
area = np.cos(np.radians(lat))[:, None] * np.ones(len(lon)) * (2 * np.pi * 6.371e6 ** 2 / len(lon)) * np.radians(np.gradient(lat))[:, None]
for name, m in [('Antarctica', lat[:, None] < -60), ('Greenland', (lat[:, None] > 59) & (lon[None, :] > 285) & (lon[None, :] < 350))]:
    print(name, 'ice-sheet area %.2f Mkm2' % ((frac * area * m).sum() / 1e12))
