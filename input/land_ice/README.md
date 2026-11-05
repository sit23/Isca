# Ice-sheet fraction on the T42 grid

`era5_land_ice_t42.nc`, field `land_ice`: the fraction of the land in each T42 cell covered by
ice sheet, for `mixed_layer_nml` `land_ice_file_name = 'era5_land_ice_t42'`. Made by
`make_land_ice_t42.py` from ERA5 monthly snow depth, 1991-2020, at 1 degree (Hersbach et al.
2020, Copernicus Climate Change Service): ERA5 keeps snow depth on glaciers and ice sheets at
10 m water equivalent, so a 1-degree land cell whose climatological snow depth is at least
1 m w.e. in every month is taken as ice sheet. Land at 1 degree is where ERA5 has no SST; the
grid and land mask are those of `exp/test_cases/simple_clouds/input/era_land_t42_filtered.nc`.
Ice-sheet area: Antarctica 13.7 Mkm2, Greenland box (59N+, 285-350E, includes some Canadian
Arctic ice caps) 2.2 Mkm2.
