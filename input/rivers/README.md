# River network

`ddm30_flowdir_cru_neva.nc`: DDM30, the global 30-minute drainage direction map of Döll and
Lehner (2002), in the ISIMIP3 version adapted to the CRU land mask (Caspian Sea excluded).
Values: 1 E, 2 SE, 3 S, 4 SW, 5 W, 6 NW, 7 N, 8 NE, 0 outlet or sink.

Source: ISIMIP repository, doi:10.48364/ISIMIP.865475 (CC BY 4.0).
Reference: Döll, P. and Lehner, B. (2002), Validation of a new global 30-min drainage
direction map, J. Hydrol. 258, 214-231.

Used by `isca.coupling.rivers` to route land runoff to real river mouths in Earth-like
coupled runs.
