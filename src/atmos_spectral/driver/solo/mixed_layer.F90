!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mixed_layer_mod

!
! Implementation of mixed layer boundary condition
!

use            fms_mod, only: set_domain, write_version_number, &
                              mpp_pe, mpp_root_pe, error_mesg, NOTE, FATAL, WARNING

use            fms_mod, only: stdlog, check_nml_error, close_file,&
                              open_namelist_file, stdout, file_exist, &
                              read_data, write_data, open_file, &
                              nullify_domain, lowercase, field_size

use  field_manager_mod, only: MODEL_ATMOS

use tracer_manager_mod, only: get_tracer_names, get_number_tracers

use      constants_mod, only: HLV, PI, RHO_CP, CP_AIR, KELVIN, & 
                              HLF_WATER, DENS_ICE, TFREEZE_WATER ! sea ice: water constants (HLF, TFREEZE are run-time variables)

use   diag_manager_mod, only: register_diag_field, register_static_field, send_data

use   time_manager_mod, only: time_type

#ifdef COLUMN_MODEL
use    column_grid_mod, only: get_deg_lon, get_deg_lat, area_weighted_global_mean
use         column_mod, only: get_surf_geopotential
use       spec_mpp_mod, only: grid_domain
#else          
use     transforms_mod, only: get_deg_lat, get_deg_lon, grid_domain, area_weighted_global_mean
! mj know about surface topography
use spectral_dynamics_mod,only: get_surf_geopotential
#endif

use      vert_diff_mod, only: surf_diff_type

use         mpp_domains_mod, only: mpp_get_global_domain !s added to enable qflux reading


! mj read SSTs
use interpolator_mod, only: interpolate_type,interpolator_init&
     &,CONSTANT,interpolator
!mj q-flux
use qflux_mod, only: qflux_init,qflux,warmpool
! coupling to an external ocean model (MOM6, ML ocean) through OASIS3-MCT
use ocean_coupler_mod, only: ocean_coupler_init, ocean_coupler_put, ocean_coupler_get, &
                             ocean_coupler_has_ice, ocean_coupler_get_ice, ocean_coupler_put_ice

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.F90,v 1.1.2.1 2013/01/24 20:35:37 pjp Exp $'

character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, mixed_layer, mixed_layer_end, albedo_calc, mixed_layer_uncoupled_pme
public :: mixed_layer_tile_albedo, mixed_layer_tile_diagnostics, mixed_layer_tile_coupling

!=================================================================================================================================

logical :: evaporation = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
logical :: load_qflux = .false.
logical :: time_varying_qflux = .false.
real    :: tconst = 305.0
real    :: delta_T = 40.0
logical :: prescribe_initial_dist = .false.
real    :: albedo_value    = 0.06

!s Surface heat capacity options
real    :: depth           = 40.0,         & !s 2013 implementation
           land_depth      = -1.,          & !mj 2013
           trop_depth      = -1.,          & !mj 2013
           trop_cap_limit  = 15.,          & !mj
           heat_cap_limit  = 60.,          & !mj
           np_cap_factor   =  1.,          & !mj
           albedo_exp      = 2.,           & !mj
           albedo_cntr     = 45.,          & !mj
           albedo_wdth     = 10.,          & !mj
           higher_albedo   = 0.10,         & !mj
           lat_glacier     = 60.,          & !mj
           land_h_capacity_prefactor = 1.0 !s where(land) heat_capcity = land_h_capacity_prefactor * depth * rho_cp

!s Surface albedo options
real    :: land_albedo_prefactor = 1.0 !s where(land) albedo = land_albedo_prefactor * albedo_value

!s Begin mj extra options
integer :: albedo_choice    = 1 ! 1->constant or following 'where(land)', 2->NH or SH step, 3->N-S symmetric step, 4->profile with albedo_exp, 5->tanh with albedo_cntr,albedo_wdth
logical :: do_qflux         = .false. !mj
logical :: do_warmpool      = .false. !mj
logical :: do_read_sst      = .false. !mj
logical :: do_sc_sst        = .false. !mj use specified SSTs
logical :: do_ape_sst       = .false. ! use the AquaPlanet Experiement (APE) sst profile.

! 'none' : surface temperature from this module's slab ocean / prescribed SSTs as usual.
! 'oasis': over points covered by an external ocean model, t_surf is the SST
!          received through OASIS3-MCT, and the surface fluxes are sent back.
!          Points not covered by the ocean (land, polar caps) keep the slab physics.
character(len=16) :: ocean_coupling = 'none'
real :: coupled_ocean_frac_threshold = 0.5 ! cells with ocean fraction above this are coupled
logical :: specify_sst_over_ocean_only = .false.
logical :: do_calc_eff_heat_cap = .true. ! assumes specified SST are off the default.

character(len=256) :: sst_file
character(len=256) :: land_option = 'none'
real,dimension(10) :: slandlon=0,slandlat=0,elandlon=-1,elandlat=-1
!s End mj extra options

logical :: specify_constant_sst = .false.
real    :: sst_prescribed_constant = 273.

character(len=256) :: qflux_file_name  = 'ocean_qflux'
character(len=256) :: qflux_field_name = 'ocean_qflux' !Only used when using a non-time-varying q-flux. Otherwise code assumes field_name = file_name. 

character(len=256) :: ice_file_name  = 'siconc_clim_amip'
real    :: ice_albedo_value = 0.7
real    :: ice_concentration_threshold = 0.5
logical :: update_albedo_from_ice = .false.
character(len=256) :: ice_albedo_method = 'step_function'

logical :: add_latent_heat_flux_anom = .false.
character(len=256) :: flux_lhe_anom_file_name  = 'INPUT/flux_lhe_anom.nc'
character(len=256) :: flux_lhe_anom_field_name = 'flux_lhe_anom'

! THERMO ICE OPTIONS - NTL 01/23
logical :: do_explicit = .false. ! explicit or implicit timestepping for ocean covered surface 
logical :: do_thermo_ice = .false. ! do thermodynamic sea ice? 
real :: thermo_ice_albedo = 0.4 ! from Feldl & Merlis - original XZ is 0.8, M. England 0.6 
real :: L_thermo_ice = DENS_ICE * HLF_WATER ! latent heat of fusion (J/m^3) 
real, parameter  :: k_thermo_ice = 2 ! conductivity of ice (W/m/K)
real, parameter  :: thermo_ice_basal_flux_const = 120 ! linear coefficient for heat flux from ocean to ice base (W/m^2/K)
real :: t_thermo_ice_base = TFREEZE_WATER ! temperature at base of ice, taken to be freshwater freezing point
real :: t_surf_freeze = TFREEZE_WATER
real :: thermo_ice_nudge_timescale = 7. * 86400. 
logical :: do_nudge_ice = .false. 
logical :: do_nudging_correction = .false.
character(len=256) :: nudge_ice_file_name = 'nudge_ice'
logical :: time_varying_nudge_ice = .false.
logical :: do_prescribe_albedo = .false.
character(len=256) :: albedo_file_name = 'albedo'
logical :: time_varying_albedo = .false.
logical :: prescribe_nh_albedo_only = .false.
logical :: read_const_correct = .true. ! set to false if doing nudging, and want to restart from simulation with no nudging
logical :: read_nudge_out = .true. 
logical :: albedo_from_nudge_file = .false.

namelist/mixed_layer_nml/ evaporation, depth, qflux_amp, qflux_width, tconst,&
                              delta_T, prescribe_initial_dist,albedo_value,  &
                              land_depth,trop_depth,                         &  !mj
                              trop_cap_limit, heat_cap_limit, np_cap_factor, &  !mj
                              do_qflux,do_warmpool,                          &  !mj
                              albedo_choice,higher_albedo,albedo_exp,        &  !mj
                              albedo_cntr,albedo_wdth,lat_glacier,           &  !mj
                              do_read_sst,do_sc_sst,sst_file,                &  !mj
                              land_option,slandlon,slandlat,                 &  !mj
                              elandlon,elandlat,                             &  !mj
                              land_h_capacity_prefactor,                     &  !s
                              land_albedo_prefactor,                         &  !s
                              load_qflux,qflux_file_name,time_varying_qflux, &
                              update_albedo_from_ice, ice_file_name,         &
                              ice_albedo_value, specify_sst_over_ocean_only, &
                              ice_concentration_threshold, ice_albedo_method,&
                              add_latent_heat_flux_anom,flux_lhe_anom_file_name,&
                              flux_lhe_anom_field_name, specify_constant_sst,&
                              sst_prescribed_constant,                       &
                              do_ape_sst, qflux_field_name,                  &
                              ocean_coupling, coupled_ocean_frac_threshold,  &
                              do_explicit, do_thermo_ice, thermo_ice_albedo, t_surf_freeze, t_thermo_ice_base, & ! NTL ice 
                              do_nudge_ice, nudge_ice_file_name, thermo_ice_nudge_timescale, time_varying_nudge_ice, &
                              do_nudging_correction, do_prescribe_albedo, time_varying_albedo, albedo_file_name, &
                              prescribe_nh_albedo_only, read_const_correct, read_nudge_out, albedo_from_nudge_file

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter, nhum
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux
     id_flux_t,            &   ! sensible heat flux at surface
     id_heat_cap,          &   ! heat capacity
     id_albedo,            &   ! mj albedo
     id_ice_conc,          &   ! st ice concentration
     id_delta_t_surf,      &
     ! NTL 01/23 thermodynamic sea ice 
     id_h_thermo_ice,      &   ! sea ice thickness
     id_t_ml,              &   ! mixed layer temperature
     id_flux_thermo_ice        ! conductive heat flux through ice

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux
     ice_concentration,     &   ! ice_concentration
     rad_lat_2d,            &   ! latitude in radians
     flux_lhe_anom, flux_q_total

real, allocatable, dimension(:)   :: deg_lat, deg_lon

real, allocatable, dimension(:,:)   ::                                        &
     gamma_t,               &   ! Used to calculate the implicit
     gamma_q,               &   ! correction to the diffusion in
     fn_t,                  &   ! the lowest layer
     fn_q,                  &   !
     en_t,                  &   !
     en_q,                  &   !
     alpha_t,               &   !
     alpha_q,               &   !
     alpha_lw,              &   !
     beta_t,                &   !
     beta_q,                &   !
     beta_lw,               &   !
     t_surf_dependence,     &   !
     corrected_flux,        &   !
     eff_heat_capacity,     &   ! Effective heat capacity
     delta_t_surf,          &   ! Increment in surface temperature
     zsurf,                 &   ! mj know about topography
     land_sea_heat_capacity,&
     sst_new,               &   ! mj input SST
     albedo_initial,        &
     ! NTL 01/23 thermodynamic sea ice 
     dFdt_surf,             &   ! d(corrected_flux)/d(t_surf), for calculation of ice sfc temperature
     delta_t_ml,            &   ! Increment in mixed layer temperature
     delta_h_thermo_ice,           &   ! Increment in sea ice thickness
     delta_t_thermo_ice,           &   ! Increment in ice (steady-state) surface temperature
     flux_thermo_ice,              &   ! Conductive heat flux through ice
     corrected_flux_noq
real, allocatable, dimension(:,:) :: h_thermo_ice_in, nudge_flux, nudge_correction, albedo_in
real :: const_nudge_correction 
type(interpolate_type),save :: nudge_ice_interp
type(interpolate_type),save :: albedo_interp

logical, allocatable, dimension(:,:) ::      land_mask

! Fractional surface tiles (idealized_moist_phys_nml surface_tiles): every cell has an
! ocean tile (1, the coupled ocean when there is one; open water when the ocean component
! has a sea-ice model), a land tile (2), a slab-ocean tile (3, ocean the coupled ocean
! model does not cover) and a sea-ice tile (4, the ice of the ocean component's sea-ice
! model, e.g. SIS2), and mixed_layer is called once for each. These hold
! the land mask the model was started with, and the heat capacity and initial albedo of
! each tile; land_mask, land_sea_heat_capacity and albedo_initial are set from them for
! the tile being stepped.
logical :: surface_tiles = .false.
logical, allocatable, dimension(:,:) :: land_mask_grid
real, allocatable, dimension(:,:,:) :: heat_cap_tile, albedo_initial_tile
real, allocatable, dimension(:,:,:) :: tile_evap   ! coupled: each tile's evaporation [kg/m2/s]
! with an ocean component that has a sea-ice model: its ice fraction (of the ocean), surface
! temperature and albedo as last received, and the albedo the radiation last used for the ice
real, allocatable, dimension(:,:) :: ice_frac_cpl, ice_ts_cpl, ice_alb_cpl, ice_alb_rad
logical :: coupled_ice = .false.

logical :: coupled = .false.
real, allocatable, dimension(:,:) :: sst_coupled, ocean_frac_coupled
real, allocatable, dimension(:,:) :: uncoupled_pme   ! P-E [kg/m2/s] at uncoupled points (see mixed_layer_uncoupled_pme)
logical, allocatable, dimension(:,:) :: coupled_mask
! with do_thermo_ice: frazil heat the ocean reports [W/m2 of Isca cell, from all the
! ocean under the cell], and the heat going into the ocean under the ice [W/m2 of
! ocean, positive down]
real, allocatable, dimension(:,:) :: frazil_coupled, ice_ocean_heat
logical, allocatable, dimension(:,:) :: ice_covered_coupled
integer :: id_ocean_frac, id_cpl_heat, id_cpl_evap, id_cpl_precip, id_cpl_frazil, id_ice_ocean_heat

!mj read sst from input file
  type(interpolate_type),save :: sst_interp
  type(interpolate_type),save :: qflux_interp
  type(interpolate_type),save :: ice_interp
  type(interpolate_type),save :: flux_lhe_anom_interp  

real inv_cp_air

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, &
                            h_thermo_ice, t_ml, const_correct, nudge_out, & ! NTL 01/23 thermodynamic sea ice 
                            bucket_depth, axes, Time, albedo, rad_lonb_2d,rad_latb_2d, land, restart_file_bucket_depth, &
                            tiles, land_full)

type(time_type), intent(in)       :: Time
real, intent(out), dimension(:,:) :: t_surf, albedo 
real, intent(out), dimension(:,:) :: h_thermo_ice, t_ml, const_correct, nudge_out ! NTL 01/23 thermodynamic sea ice 
real, intent(out), dimension(:,:,:) :: bucket_depth
integer, intent(in), dimension(4) :: axes
real, intent(in), dimension(:,:) :: rad_lonb_2d, rad_latb_2d
integer, intent(in) :: is, ie, js, je, num_levels

logical, intent(in), dimension(:,:) :: land
logical, intent(in)                 :: restart_file_bucket_depth
logical, intent(in), optional       :: tiles   ! fractional surface tiles (see surface_tiles)
logical, intent(in), dimension(:,:), optional :: land_full   ! with tiles: cells with no ocean tile

integer :: j
real    :: rad_qwidth
integer:: ierr, io, unit, num_tr, n, global_num_lon, global_num_lat
character(32) :: tr_name

! mj shallower ocean in tropics, land-sea contrast
 real :: trop_capacity,land_capacity,lon,lat,loc_cap
 integer :: i,k
integer, dimension(4) :: siz
character(len=12) :: ctmp1='     by     ', ctmp2='     by     '

if(module_is_initialized) return

call write_version_number(version, tagname)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=mixed_layer_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'mixed_layer_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=mixed_layer_nml)

call get_number_tracers (MODEL_ATMOS, num_prog=num_tr)
do n = 1,num_tr
   call get_tracer_names( MODEL_ATMOS, n, tr_name )
   if(lowercase(tr_name)=='sphum') then
      nhum = n
   endif
enddo

allocate(rad_lat_2d              (is:ie, js:je))
allocate(ocean_qflux             (is:ie, js:je))
allocate(ice_concentration       (is:ie, js:je))
allocate(flux_lhe_anom           (is:ie, js:je))
allocate(flux_q_total            (is:ie, js:je))
allocate(deg_lat                 (js:je))
allocate(deg_lon                 (is:ie))
allocate(gamma_t                 (is:ie, js:je))
allocate(gamma_q                 (is:ie, js:je))
allocate(en_t                    (is:ie, js:je))
allocate(en_q                    (is:ie, js:je))
allocate(fn_t                    (is:ie, js:je))
allocate(fn_q                    (is:ie, js:je))
allocate(alpha_t                 (is:ie, js:je))
allocate(alpha_q                 (is:ie, js:je))
allocate(alpha_lw                (is:ie, js:je))
allocate(beta_t                  (is:ie, js:je))
allocate(beta_q                  (is:ie, js:je))
allocate(beta_lw                 (is:ie, js:je))
allocate(delta_t_surf            (is:ie, js:je))
allocate(eff_heat_capacity       (is:ie, js:je))
allocate(corrected_flux          (is:ie, js:je))
allocate(corrected_flux_noq      (is:ie, js:je)) ! NTL 01/23 thermodynamic sea ice 
allocate(t_surf_dependence       (is:ie, js:je))
allocate (albedo_initial         (is:ie, js:je))
allocate(land_sea_heat_capacity  (is:ie, js:je))
allocate(zsurf                   (is:ie, js:je))
allocate(sst_new                 (is:ie, js:je))
allocate(land_mask                 (is:ie, js:je)); land_mask=land
allocate(land_mask_grid            (is:ie, js:je)); land_mask_grid=land
if (present(tiles)) surface_tiles = tiles
if (surface_tiles) then
  allocate(heat_cap_tile(is:ie, js:je, 4), albedo_initial_tile(is:ie, js:je, 4))
  heat_cap_tile = depth*RHO_CP ; albedo_initial_tile = albedo_value
  if (update_albedo_from_ice .or. do_nudge_ice .or. do_nudging_correction .or. add_latent_heat_flux_anom) &
    call error_mesg('mixed_layer', 'surface_tiles cannot be used with update_albedo_from_ice, '// &
                    'do_nudge_ice, do_nudging_correction or add_latent_heat_flux_anom', FATAL)
  if (trim(land_option) /= 'input' .and. trim(land_option) /= 'none') &
    call error_mesg('mixed_layer', 'surface_tiles needs land_option = "input" or "none"', FATAL)
endif

! NTL 01/23 thermodynamic sea ice 
allocate(dFdt_surf               (is:ie, js:je))
allocate(delta_t_ml              (is:ie, js:je))
allocate(delta_h_thermo_ice      (is:ie, js:je))
allocate(delta_t_thermo_ice      (is:ie, js:je))
allocate(flux_thermo_ice         (is:ie, js:je))
allocate(h_thermo_ice_in         (is:ie, js:je))
allocate(albedo_in               (is:ie, js:je))
allocate(nudge_flux              (is:ie, js:je))
allocate(nudge_correction        (is:ie, js:je))
!
!see if restart file exists for the surface temperature
!


! latitude will be needed for oceanic q flux
!s Moved up slightly so that rad_lat_2d can be used in initial temperature distribution if necessary.

if (do_thermo_ice .and. (.not. do_explicit)) then 
  call error_mesg('mixed_layer','setting do_explicit .true. for consistency with sea ice', WARNING)
  do_explicit = .true.
endif 


call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo
!mj get lon for land_sea_heat_capacity
call get_deg_lon(deg_lon)

!s Adding MiMA options
   ! if(do_sc_sst) do_read_sst = .true.
   trop_capacity   = trop_depth*RHO_CP
   land_capacity   = land_depth*RHO_CP
   if(trop_capacity .le. 0.) trop_capacity = depth*RHO_CP
   if(land_capacity .le. 0.) land_capacity = depth*RHO_CP
!s End MiMA options

    !mj read fixed SSTs
    if( do_read_sst ) then
        call interpolator_init( sst_interp, trim(sst_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    endif




if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
   if (do_thermo_ice) then 
       call read_data(trim('INPUT/mixed_layer.res'), 'h_thermo_ice', h_thermo_ice, grid_domain)
       call read_data(trim('INPUT/mixed_layer.res'), 't_ml', t_ml, grid_domain)
       if (read_const_correct) then 
         call read_data(trim('INPUT/mixed_layer.res'), 'const_correct', const_correct, grid_domain)
       else 
         const_correct(:,:) = 0.0 
       endif 
       if (read_nudge_out) then 
         call read_data(trim('INPUT/mixed_layer.res'), 'nudge_out', nudge_out, grid_domain)
       else 
         nudge_out(:,:) = 0.0 
       endif 
       call read_data(trim('INPUT/mixed_layer.res'), 'albedo', albedo, grid_domain)
   else if (do_prescribe_albedo) then 
       call read_data(trim('INPUT/mixed_layer.res'), 'albedo', albedo, grid_domain)
   endif
   
   if (restart_file_bucket_depth) then
       call read_data(trim('INPUT/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)
   endif

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)

else if( do_read_sst ) then !s Added so that if we are reading sst values then we can restart using them.

   call interpolator( sst_interp, Time, t_surf, trim(sst_file) )

elseif (prescribe_initial_dist) then
!  call error_mesg('mixed_layer','mixed_layer restart file not found - initializing from prescribed distribution', WARNING)

    
    t_surf(:,:) = tconst - delta_T*((3.*sin(rad_lat_2d)**2.)-1.)/3.
    if (do_thermo_ice) then  ! NTL 01/23 thermodynamic sea ice 
        h_thermo_ice(:,:) = 0.0 
        t_ml(:,:) = t_surf(:,:)
        const_correct(:,:) = 0.0
        nudge_out(:,:) = 0.0
        albedo(:,:) = albedo_value
        where(land) albedo(:,:) = land_albedo_prefactor * albedo_value
        albedo_initial(:,:) = albedo (:,:)
    else if (do_prescribe_albedo) then 
        albedo(:,:) = albedo_value
        where(land) albedo(:,:) = land_albedo_prefactor * albedo_value
        albedo_initial(:,:) = albedo (:,:)
    endif
else

  call error_mesg('mixed_layer','mixed_layer restart file not found - initializing from lowest model level temp', WARNING)

endif


if(trim(ice_albedo_method) == 'ramp_function') then
  call error_mesg('mixed_layer','Alternative method ramp_function used for ice albedo output.', NOTE)
endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')
id_heat_cap = register_static_field(mod_name, 'ml_heat_cap',        &
                                 axes(1:2), 'mixed layer heat capacity','joules/m^2/deg C')
id_delta_t_surf = register_diag_field(mod_name, 'delta_t_surf',        &
                                 axes(1:2), Time, 'change in sst','K')
if (update_albedo_from_ice .or. do_thermo_ice .or. do_prescribe_albedo) then
    id_albedo = register_diag_field(mod_name, 'albedo',    &
                                 axes(1:2), Time, 'surface albedo', 'none')
    if (do_thermo_ice) then 
    ! NTL 01/23 thermodynamic sea ice 
        id_h_thermo_ice = register_diag_field(mod_name, 'h_thermo_ice',        &
                                        axes(1:2), Time, 'sea ice thickness','m')
        id_t_ml = register_diag_field(mod_name, 't_ml',        &
                                        axes(1:2), Time, 'mixed layer tempeature','K')
        id_flux_thermo_ice = register_diag_field(mod_name, 'flux_thermo_ice',        &
                                        axes(1:2), Time, 'conductive heat flux through sea ice','watts/m2')
    else                                   
        id_ice_conc = register_diag_field(mod_name, 'ice_conc',    &
                                     axes(1:2), Time, 'ice_concentration', 'none')
    endif 
else
    id_albedo = register_static_field(mod_name, 'albedo',    &
                                 axes(1:2), 'surface albedo', 'none')
endif

ocean_qflux = 0.

! load Q flux
if (load_qflux) then

    if (time_varying_qflux) then
       call interpolator_init( qflux_interp, trim(qflux_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    else

       if(file_exist(trim('INPUT/'//qflux_file_name//'.nc'))) then
         call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
         call field_size(trim('INPUT/'//qflux_file_name//'.nc'), trim(qflux_field_name), siz)
         if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
           call read_data(trim('INPUT/'//qflux_file_name//'.nc'), trim(qflux_field_name), ocean_qflux, grid_domain)
         else
           write(ctmp1(1: 4),'(i4)') siz(1)
           write(ctmp1(9:12),'(i4)') siz(2)
           write(ctmp2(1: 4),'(i4)') global_num_lon
           write(ctmp2(9:12),'(i4)') global_num_lat
           call error_mesg ('get_qflux','Qflux file contains data on a '// &
                  ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
         endif
       else
         call error_mesg('get_qflux','load_qflux="'//trim('True')//'"'// &
                         ' but '//trim(qflux_file_name)//' does not exist', FATAL)
       endif

    endif
endif


h_thermo_ice_in = 999. ! this is used to determine whether nudging should occur 
                       ! setting a large number ensures this will do nothing if not initialised from external file
! load ice
if (do_nudge_ice .or. albedo_from_nudge_file) then
    if (time_varying_nudge_ice .or. albedo_from_nudge_file) then
       call interpolator_init( nudge_ice_interp, trim(nudge_ice_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    else

       if(file_exist(trim('INPUT/'//nudge_ice_file_name//'.nc'))) then
         call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
         call field_size(trim('INPUT/'//nudge_ice_file_name//'.nc'), trim(nudge_ice_file_name), siz)
         if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
           call read_data(trim('INPUT/'//nudge_ice_file_name//'.nc'), trim(nudge_ice_file_name), h_thermo_ice_in, grid_domain)
         else
           write(ctmp1(1: 4),'(i4)') siz(1)
           write(ctmp1(9:12),'(i4)') siz(2)
           write(ctmp2(1: 4),'(i4)') global_num_lon
           write(ctmp2(9:12),'(i4)') global_num_lat
           call error_mesg ('get_ice','nudge_ice file contains data on a '// &
                  ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
         endif
       else
         call error_mesg('get_ice','do_nudge_ice="'//trim('True')//'"'// &
                         ' but '//trim(nudge_ice_file_name)//' does not exist', FATAL)
       endif

    endif
endif

albedo_in = 0.0 
if (do_prescribe_albedo) then  
  if (time_varying_albedo) then
       call interpolator_init( albedo_interp, trim(albedo_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    else

       if(file_exist(trim('INPUT/'//albedo_file_name//'.nc'))) then
         call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
         call field_size(trim('INPUT/'//albedo_file_name//'.nc'), trim(albedo_file_name), siz)
         if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
           call read_data(trim('INPUT/'//albedo_file_name//'.nc'), trim(albedo_file_name), albedo_in, grid_domain)
         else
           write(ctmp1(1: 4),'(i4)') siz(1)
           write(ctmp1(9:12),'(i4)') siz(2)
           write(ctmp2(1: 4),'(i4)') global_num_lon
           write(ctmp2(9:12),'(i4)') global_num_lat
           call error_mesg ('get_albedo','albedo file contains data on a '// &
                  ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
         endif
       else
         call error_mesg('get_albedo','do_prescribe_albedo="'//trim('True')//'"'// &
                         ' but '//trim(albedo_file_name)//' does not exist', FATAL)
       endif

    endif
endif


!s Adding MiMA options for qfluxes.

if ( do_qflux .or. do_warmpool) then
   call qflux_init
!mj q-flux as in Merlis et al (2013) [Part II]
   if ( do_qflux ) call qflux(rad_latb_2d(1,:),ocean_qflux)
!mj q-flux to create a tropical temperature perturbation
   if ( do_warmpool) call warmpool(rad_lonb_2d(:,1),rad_latb_2d(1,:),ocean_qflux)
endif

!s End MiMA options for qfluxes

if (add_latent_heat_flux_anom) then
       call interpolator_init( flux_lhe_anom_interp, trim(flux_lhe_anom_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
endif



inv_cp_air = 1.0 / CP_AIR

!s Prescribe albedo distribution here so that it will be the same in both two_stream_gray later and rrtmg radiation.

if ((.not. do_thermo_ice) .and. (.not. do_prescribe_albedo)) then ! NTL 01/23 thermodynamic sea ice -- this is set / loaded earlier 

    albedo(:,:) = albedo_value

    if(trim(land_option) .eq. 'input') then

    where(land) albedo = land_albedo_prefactor * albedo

    endif

    !mj MiMA albedo choices.
    select case (albedo_choice)
      case (2) ! higher_albedo northward (lat_glacier>0) or southward (lat_glacier <0 ) of lat_glacier
        do j = 1, size(t_surf,2)
          lat = deg_lat(js+j-1)
          ! mj SH or NH only
          if ( lat_glacier .ge. 0. ) then
             if ( lat > lat_glacier ) then
                albedo(:,j) = higher_albedo
             else
                albedo(:,j) = albedo_value
             endif
          else
             if ( lat < lat_glacier ) then
                albedo(:,j) = higher_albedo
             else
                albedo(:,j) = albedo_value
             endif
          endif
        enddo
      case (3) ! higher_albedo poleward of lat_glacier
        do j = 1, size(t_surf,2)
          lat = deg_lat(js+j-1)
          if ( abs(lat) > lat_glacier ) then
            albedo(:,j) = higher_albedo
          else
            albedo(:,j) = albedo_value
          endif
        enddo
      case (4) ! exponential increase with albedo_exp
         do j = 1, size(t_surf,2)
           lat = abs(deg_lat(js+j-1))
           albedo(:,j) = albedo_value + (higher_albedo-albedo_value)*(lat/90.)**albedo_exp
         enddo
      case (5) ! tanh increase around albedo_cntr with albedo_wdth
         do j = 1, size(t_surf,2)
           lat = abs(deg_lat(js+j-1))
           albedo(:,j) = albedo_value + (higher_albedo-albedo_value)*&
                 0.5*(1+tanh((lat-albedo_cntr)/albedo_wdth))
         enddo
    end select

    albedo_initial=albedo
    if (surface_tiles) then
      ! as above, with the whole tile ocean or land: only albedo_choice 1 tells them apart
      albedo_initial_tile(:,:,1) = albedo
      albedo_initial_tile(:,:,2) = albedo
      if (albedo_choice == 1) then
        albedo_initial_tile(:,:,1) = albedo_value
        if (trim(land_option) .eq. 'input') albedo_initial_tile(:,:,2) = land_albedo_prefactor * albedo_value
      endif
      albedo_initial_tile(:,:,3) = albedo_initial_tile(:,:,1)
      albedo_initial_tile(:,:,4) = albedo_initial_tile(:,:,1)
    endif

    if (update_albedo_from_ice) then
        call interpolator_init( ice_interp, trim(ice_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
        call read_ice_conc(Time)
        call albedo_calc(albedo,Time)
    else
        if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo )
    endif
else ! NTL 01/23 thermodynamic sea ice 
    if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo, Time) 
endif ! if(do_thermo_ice) 

! Note: do_calc_eff_heat_cap is true by default and control when the surface 
! heat capacity is calculated (land and ocean). 
if (do_sc_sst) then
    ! if using specified sst then do not calc the heat capacity
    if (specify_sst_over_ocean_only) then
        ! unless the heat capacity over land is needed
        do_calc_eff_heat_cap = .true.
    else
        do_calc_eff_heat_cap = .false.
    endif
else
    if (do_ape_sst) then
        ! if using specified sst without land, do not calc the heat capacity
        do_calc_eff_heat_cap = .false.
    end if
endif



!s begin surface heat capacity calculation
if (do_calc_eff_heat_cap) then
    land_sea_heat_capacity = depth*RHO_CP
    if(trim(land_option) .ne. 'input') then
         if ( trop_capacity .ne. depth*RHO_CP .or. np_cap_factor .ne. 1. ) then !s Lines above make trop_capacity=depth*RHO_CP if trop_capacity set to be < 0.
            do j=js,je
               lat = deg_lat(j)
               if ( lat .gt. 0. ) then
                  loc_cap = depth*RHO_CP*np_cap_factor
               else
                  loc_cap = depth*RHO_CP
               endif
               if ( abs(lat) .lt. trop_cap_limit ) then
                  land_sea_heat_capacity(:,j) = trop_capacity
               elseif ( abs(lat) .lt. heat_cap_limit ) then
                  land_sea_heat_capacity(:,j) = trop_capacity*(1.-(abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)) + (abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)*loc_cap
               elseif ( lat .gt. heat_cap_limit ) then
                  land_sea_heat_capacity(:,j) = loc_cap
               end if
            enddo
         endif
! mj land heat capacity function of surface topography
         if(trim(land_option) .eq. 'zsurf')then
            call get_surf_geopotential(zsurf)
            where ( zsurf .gt. 10. ) land_sea_heat_capacity = land_capacity
         endif
! mj land heat capacity given through ?landlon, ?landlat
         if(trim(land_option) .eq. 'lonlat')then
            do j=js,je
           lat = deg_lat(j)
               do i=is,ie
                  lon = deg_lon(i)
                  do k=1,size(slandlat)
                     if ( lon .ge. slandlon(k) .and. lon .le. elandlon(k) &
                          &.and. lat .ge. slandlat(k) .and. lat .le. elandlat(k) )then
                        land_sea_heat_capacity(i,j) = land_capacity
                     endif
                  enddo
               enddo
            enddo
         endif
    else  !trim(land_option) .eq. 'input'
        if (surface_tiles) then
          heat_cap_tile(:,:,1) = land_sea_heat_capacity
          heat_cap_tile(:,:,2) = land_h_capacity_prefactor*land_sea_heat_capacity
        endif
        where(land) land_sea_heat_capacity = land_h_capacity_prefactor*land_sea_heat_capacity
    endif !end of if (trim(land_option) .ne. 'input')
    if (surface_tiles .and. trim(land_option) .ne. 'input') then
      heat_cap_tile(:,:,1) = land_sea_heat_capacity
      heat_cap_tile(:,:,2) = land_sea_heat_capacity
    endif
    if (surface_tiles) heat_cap_tile(:,:,3) = heat_cap_tile(:,:,1)
    if (surface_tiles) heat_cap_tile(:,:,4) = heat_cap_tile(:,:,1)
endif !end of if(.not.do_sc_sst)

if ( id_heat_cap > 0 ) used = send_data ( id_heat_cap, land_sea_heat_capacity )
!s end surface heat capacity calculation

if (trim(ocean_coupling) == 'oasis') then
  coupled = .true.
  allocate(sst_coupled(is:ie, js:je), ocean_frac_coupled(is:ie, js:je), coupled_mask(is:ie, js:je))
  allocate(uncoupled_pme(is:ie, js:je)); uncoupled_pme = 0.
  sst_coupled = t_surf
  ocean_frac_coupled = 0.
  coupled_mask = .false.
  call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
  ! with tiles every cell with an ocean tile receives the ocean's fields
  if (surface_tiles .and. present(land_full)) then
    call ocean_coupler_init(Time, is, ie, js, je, global_num_lon, global_num_lat, &
                            rad_lonb_2d, rad_latb_2d, land_full)
  else
    call ocean_coupler_init(Time, is, ie, js, je, global_num_lon, global_num_lat, &
                            rad_lonb_2d, rad_latb_2d, land)
  endif
  id_ocean_frac = register_diag_field(mod_name, 'coupled_ocean_frac', axes(1:2), Time, &
                                      'fraction of cell covered by the coupled ocean', 'none')
  ! exactly what is sent to the ocean (all points; the ocean uses the coupled ones)
  id_cpl_heat = register_diag_field(mod_name, 'cpl_net_heat', axes(1:2), Time, &
                  'net surface heat flux sent to the ocean (SW+LW-SH-L*E), positive down', 'W/m2')
  id_cpl_evap = register_diag_field(mod_name, 'cpl_evap', axes(1:2), Time, &
                  'evaporation sent to the ocean', 'kg/m2/s')
  id_cpl_precip = register_diag_field(mod_name, 'cpl_precip', axes(1:2), Time, &
                  'precipitation sent to the ocean', 'kg/m2/s')
  allocate(frazil_coupled(is:ie, js:je), ice_ocean_heat(is:ie, js:je), ice_covered_coupled(is:ie, js:je))
  frazil_coupled = 0. ; ice_ocean_heat = 0. ; ice_covered_coupled = .false.
  if (surface_tiles) then
    allocate(tile_evap(is:ie, js:je, 4)); tile_evap = 0.
    coupled_ice = ocean_coupler_has_ice()
    if (coupled_ice) then
      if (do_thermo_ice) call error_mesg('mixed_layer', 'the ocean component has a sea-ice '// &
          'model: do_thermo_ice must be off', FATAL)
      allocate(ice_frac_cpl(is:ie, js:je), ice_ts_cpl(is:ie, js:je), ice_alb_cpl(is:ie, js:je), &
               ice_alb_rad(is:ie, js:je))
      ice_frac_cpl = 0. ; ice_ts_cpl = t_surf ; ice_alb_cpl = albedo_value ; ice_alb_rad = albedo_value
    endif
  else if (ocean_coupler_has_ice()) then
    call error_mesg('mixed_layer', 'the ocean component has a sea-ice model: '// &
        'it needs surface_tiles', FATAL)
  endif
  if (do_thermo_ice) then
    if (update_albedo_from_ice) call error_mesg('mixed_layer', 'do_thermo_ice and the prescribed-ice '// &
        'lid (update_albedo_from_ice) cannot both be used with a coupled ocean', FATAL)
    id_cpl_frazil = register_diag_field(mod_name, 'cpl_frazil', axes(1:2), Time, &
                    'frazil heat received from the ocean (turned into sea ice)', 'W/m2')
    id_ice_ocean_heat = register_diag_field(mod_name, 'ice_ocean_heat', axes(1:2), Time, &
                    'heat flux into the coupled ocean under sea ice (melt minus basal flux)', 'W/m2 of ocean')
  endif
else if (trim(ocean_coupling) /= 'none') then
  call error_mesg('mixed_layer','"'//trim(ocean_coupling)//'" is not a valid ocean_coupling. '// &
                  'Choices are: none or oasis.', FATAL)
endif

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine mixed_layer (                                               &
     Time,                                                             &
     Time_next,                                                        &
     js, je,                                                           &
     t_surf,                                                           &
     ! NTL 01/23 thermodynamic sea ice
     h_thermo_ice,                                                     &
     t_ml,                                                             &
     const_correct, &
     nudge_out, &
     ! end ice addition
     flux_t,                                                           &
     flux_q,                                                           &
     flux_r,                                                           &
     dt,                                                               &
     net_surf_sw_down,                                                 &
     surf_lw_down,                                                     &
     Tri_surf,                                                         &
     dhdt_surf,                                                        &
     dedt_surf,                                                        &
     dedq_surf,                                                        &
     drdt_surf,                                                        &
     dhdt_atm,                                                         &
     dedq_atm,                                                         &
     albedo_out,                                                       &
     flux_u, flux_v, precip, psurf, tile, dhdt_atm_grid, dedq_atm_grid, t_atm)

! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)         :: Time, Time_next
integer, intent(in)                 :: js, je
real, intent(in), dimension(:,:)    :: net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:)    :: flux_t, flux_q, flux_r
real, intent(inout), dimension(:,:) :: t_surf
real, intent(inout), dimension(:,:) :: h_thermo_ice, t_ml, const_correct, nudge_out ! NTL 01/23 thermodynamic sea ice
real, intent(in), dimension(:,:)    :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm, dedq_atm
real, intent(in)                    :: dt
real, intent(out), dimension(:,:)   :: albedo_out
type(surf_diff_type), intent(inout) :: Tri_surf
! needed only when coupled to an external ocean: surface stress on the atmosphere [Pa],
! precipitation reaching the surface [kg/m2/s] and surface pressure [Pa]
real, intent(in), dimension(:,:), optional :: flux_u, flux_v, precip, psurf
real, intent(in), dimension(:,:), optional :: t_atm   ! lowest-level air temperature (ice-tile snow)
! with surface_tiles: the tile being stepped (1 ocean, 2 land, 3 slab ocean), and the cell's total
! d(sensible heat)/d(T_atm) and d(evaporation)/d(q_atm) over all tiles, which set the
! implicit response of the lowest atmospheric level. The caller combines the tiles'
! increments to the lowest level, Tri_surf%delta_t and delta_tr, weighted by area.
integer, intent(in), optional :: tile
real, intent(in), dimension(:,:), optional :: dhdt_atm_grid, dedq_atm_grid
logical :: ocean_tile, grid_diags, ice_tile
real, dimension(size(t_surf,1),size(t_surf,2)) :: sens_out, snow_out, evap_out, lwnet_out, swnet_out, dts_flux
real, dimension(size(t_surf,1),size(t_surf,2)) :: h_ice_old, t_surf_old, frac_cpl, basal_flux, melt_excess
real, dimension(size(t_surf,1),size(t_surf,2)) :: h_old_save, t_surf_save, frazil_here
logical, dimension(size(t_surf,1),size(t_surf,2)) :: cpl
logical :: sst_received
logical, dimension(size(land_mask,1),size(land_mask,2)) :: land_ice_mask

!local variables

integer :: j

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

! Without tiles this is the whole surface, with them one tile: the ocean tile does the
! coupling, and diagnostics of the whole cell are sent by mixed_layer_tile_diagnostics.
ocean_tile = .true.
grid_diags = .true.
ice_tile = .false.
if (present(tile)) then
  if (.not. surface_tiles) call error_mesg('mixed_layer', 'tile given but surface_tiles is off', FATAL)
  ocean_tile = (tile == 1)
  ice_tile = (tile == 4)
  if (ice_tile .and. .not. coupled_ice) call error_mesg('mixed_layer', &
      'sea-ice tile without a sea-ice model in the ocean component', FATAL)
  grid_diags = .false.
  land_mask = (tile == 2)
  land_sea_heat_capacity = heat_cap_tile(:,:,tile)
  if (allocated(albedo_initial_tile)) albedo_initial = albedo_initial_tile(:,:,tile)
  if (coupled .and. .not. ocean_tile) coupled_mask = .false.
else if (surface_tiles) then
  call error_mesg('mixed_layer', 'surface_tiles is on, so mixed_layer needs tile', FATAL)
endif

if(update_albedo_from_ice) then
  call read_ice_conc(Time_next)
  land_ice_mask=.false.
  where(land_mask.or.(ice_concentration.gt.ice_concentration_threshold))
    land_ice_mask=.true.
  end where
else
  land_ice_mask=land_mask
endif

if (.not. do_thermo_ice) then ! NTL 01/23 thermodynamic sea ice
    call albedo_calc(albedo_out,Time_next)
endif 

!s Add latent heat flux anomalies before any of the calculations take place

if (add_latent_heat_flux_anom) then
    call interpolator( flux_lhe_anom_interp, Time, flux_lhe_anom, trim(flux_lhe_anom_file_name) )
    flux_q_total = flux_q + flux_lhe_anom
else
    flux_q_total = flux_q
endif


! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps

! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature

if (present(dhdt_atm_grid) .and. present(dedq_atm_grid)) then
  gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm_grid * inv_cp_air))
  gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_tr(:,:,nhum) + dedq_atm_grid))
else
  gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
  gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_tr(:,:,nhum) + dedq_atm))
endif

fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_tr(:,:,nhum) + Tri_surf%dtmass * flux_q_total)

en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf

!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q_total + dedq_atm * fn_q
alpha_lw = flux_r

beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

! If time-varying qflux then update value
if(load_qflux.and.time_varying_qflux) then
         call interpolator( qflux_interp, Time, ocean_qflux, trim(qflux_file_name) )

     if(update_albedo_from_ice) then
          where (land_ice_mask) ocean_qflux=0.
     endif

endif

if((do_nudge_ice.and.time_varying_nudge_ice).or. albedo_from_nudge_file) then 
         call interpolator( nudge_ice_interp, Time_next, h_thermo_ice_in, trim(nudge_ice_file_name) )
endif

if(do_prescribe_albedo.and.time_varying_albedo) then 
         call interpolator( albedo_interp, Time, albedo_in, trim(albedo_file_name) )
endif






!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw - ocean_qflux
t_surf_dependence = beta_t * CP_AIR + beta_lw







nudge_out = 0.0
if (do_nudge_ice) then 
  where ((h_thermo_ice .gt. 0.0) .and. (h_thermo_ice_in .le. 0)) 
    nudge_out = h_thermo_ice / thermo_ice_nudge_timescale * L_thermo_ice 
  endwhere 
  corrected_flux = corrected_flux - nudge_out 
endif 
  


! if (prescribe_nh_albedo_only) then 
!   where (rad_lat_2d .gt. 0) 

const_nudge_correction = 0.0 
const_correct = 0.0
if (do_nudging_correction) then 
  const_nudge_correction = -area_weighted_global_mean(nudge_out)
  if (.not. (const_nudge_correction .eq. 0.)) then 
    where (nudge_out .eq. 0.) 
      const_correct(:,:) = const_nudge_correction 
    endwhere 
    const_correct = const_correct * const_nudge_correction / area_weighted_global_mean(const_correct)
    corrected_flux = corrected_flux - const_correct
  endif 
endif 

if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

if (coupled .and. ocean_tile) then
  ! SST from the ocean component: only changes at coupling times, so between
  ! them delta_t_surf is zero over the coupled points.
  call ocean_coupler_get(Time, sst_coupled, ocean_frac_coupled, sst_received, frazil_coupled)
  if (coupled_ice) call ocean_coupler_get_ice(Time, ice_frac_cpl, ice_ts_cpl, ice_alb_cpl)
  coupled_mask = (ocean_frac_coupled > coupled_ocean_frac_threshold) .and. .not. land_ice_mask
  ! with tiles the ocean tile is the ocean model's ocean wherever there is any
  if (surface_tiles) coupled_mask = (ocean_frac_coupled > 0.)
  ! open water takes the ocean's SST; over sea ice (do_thermo_ice) t_surf is the ice surface
  where (coupled_mask .and. h_thermo_ice <= 0.)
    delta_t_surf = sst_coupled - t_surf
    t_surf = sst_coupled
  end where
  ! any q-flux is the ocean model's job over coupled points
  where (coupled_mask) corrected_flux = corrected_flux + ocean_qflux
  if (id_ocean_frac > 0) used = send_data(id_ocean_frac, ocean_frac_coupled, Time)
endif
if (do_thermo_ice) then 
    ! NTL 01/23 thermodynamic sea ice
    corrected_flux_noq = corrected_flux + ocean_qflux 
    flux_thermo_ice = 0 ! Initialize conductive heat flux through sea ice
    nudge_flux = 0. ! Initialize nudging flux 
    nudge_correction = 0. ! Initializing nudging correction 
    const_nudge_correction = 0.
    ! for calculation of ice sfc temperature
    dFdt_surf = dhdt_surf + drdt_surf + HLV * (dedt_surf) ! d(corrected_flux)/d(T_surf)
endif 

!s Surface heat_capacity calculation based on that in MiMA by mj

if(do_sc_sst) then !mj sst read from input file
         ! read at the new time, as that is what we are stepping to
      if (specify_constant_sst) then
         sst_new = sst_prescribed_constant
      else
         call interpolator( sst_interp, Time_next, sst_new, trim(sst_file) )
      endif

         if(specify_sst_over_ocean_only) then
            where (.not.land_ice_mask) delta_t_surf = sst_new - t_surf
            where (.not.land_ice_mask) t_surf = t_surf + delta_t_surf			 
         else
            delta_t_surf = sst_new - t_surf
            t_surf = t_surf + delta_t_surf
         endif
	     
end if

if (do_ape_sst) then 
    !
    ! AquaPlanet Experiment protocol (APE) from 
    ! Williams et al 2012 tech report: "The APE Atlas"
    !
    ! use analytic form for setting SST at each timestep.
    ! see appendix equation 1 of Neale and Hoskins 2000
    !     "A standard test for AGCMs including their 
    !      physical parametrizations: I: the proposal"
    ! using a constant longitude.
    do j=js,je   
        if ( (rad_lat_2d(1,j) .gt. -PI/3.) .and. (rad_lat_2d(1,j) .lt. PI/3.) ) then 
            ! between 60N-60S
            sst_new(:,j) = KELVIN+( 27.0*( 1. - (sin( 3./2. * rad_lat_2d(:,j) )**2 ) ))  
            !write(6,*) 'SST profile', rad_lat_2d(1,j)*180/PI, ' j:', j, ' sst:', sst_new(1,j)
        else
            ! from 60N/S to pole
            sst_new(:,j) = KELVIN
            !write(6,*) 'SST is zero', rad_lat_2d(1,j)*180/PI, ' j:', j, ' sst:', sst_new(1,j)
        endif
    enddo
    delta_t_surf = sst_new - t_surf
    t_surf = sst_new
endif


if (ice_tile) then
    ! sea ice of the ocean component's ice model: its surface temperature and albedo are
    ! prescribed here, and the ice model steps them
    where (ice_frac_cpl > 0.)
      delta_t_surf = ice_ts_cpl - t_surf
      t_surf = ice_ts_cpl
      albedo_out = ice_alb_cpl
    elsewhere
      delta_t_surf = 0.
    end where
else if (.not. do_thermo_ice) then ! NTL 01/23 thermodynamic sea ice switch 

    if (do_calc_eff_heat_cap) then
      !s use the land_sea_heat_capacity calculated in mixed_layer_init

        ! Now update the mixed layer surface temperature using an implicit step
        !
        if (do_explicit) then ! NTL 01/23 thermodynamic sea ice -- add explicit option for consistency with thermo ice 
            eff_heat_capacity = land_sea_heat_capacity
        else
            eff_heat_capacity = land_sea_heat_capacity + t_surf_dependence * dt !s need to investigate how this works
        endif

        if (any(eff_heat_capacity .eq. 0.0))  then
          write(*,*) 'mixed_layer: error', eff_heat_capacity
          call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
        end if
    if(do_sc_sst.and.specify_sst_over_ocean_only) then
        where (land_ice_mask) delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
        where (land_ice_mask) t_surf = t_surf + delta_t_surf             
    else if (coupled) then
        where (.not. coupled_mask) delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
        where (.not. coupled_mask) t_surf = t_surf + delta_t_surf
    else
        delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
        t_surf = t_surf + delta_t_surf
    endif


    endif !s end of if(do_sc_sst).
    
else if (do_thermo_ice) then 
    h_old_save = h_thermo_ice
    t_surf_save = t_surf
    where (land_ice_mask) ! name misleading, for thermo sea ice, land_ice_mask is just land mask 
        delta_t_ml = - corrected_flux * dt / land_sea_heat_capacity
        delta_h_thermo_ice  = 0
        ! t_surf and t_ml are equal (they differ only when mixed layer is ice-covered)
        t_surf = t_ml + delta_t_ml
        ! = update state =
        t_ml = t_ml + delta_t_ml
        h_thermo_ice = h_thermo_ice + delta_h_thermo_ice
    elsewhere 
    
        ! NTL 01/23 thermodynamic sea ice
        ! COPIED FROM Sea ice by Ian Eisenman, XZ 02/2018
        ! === (3) surface temperature increment is calculated with
        ! explicit forward time step using flux corrected for implicit atm;
        ! next, (4) surface state is updated  ===

        ! ======================================================================
        !
        !          Ocean mixed layer and sea ice model equations
        !
        ! ======================================================================
        !
        ! Mixed layer temperature is evolved is t_ml, and sea ice thickness is
        ! h_ice. Atmosphere cares about surface temparature (t_surf). Where
        ! h_ice=0, t_surf=t_ml, but where h_ice>0, t_surf=t_ice which is the
        ! ice surface temperature (assumed to be in steady-state: "zero layer
        ! model"). So h_ice, t_ml, and t_surf are all saved at each time step
        ! (t_ice does not need to be saved).
        !
        ! This model is equivalent to the Semtner (1976) "zero-layer model",
        ! with several differences (no snow, ice latent heat is same at base
        ! as surface, ice surface temperature calculated by including sensible
        ! and latent heat derivatives in dF/dT_surf, inclusion of frazil
        ! growth).

      ! calculate values of delta_h_ice, delta_t_ml, and delta_t_ice
         
         
         
         
         where ( h_thermo_ice .le. 0 ) ! = ice-free = [ should be equivalent to use .eq. instead of .le. ]
           delta_t_ml = - ( corrected_flux_noq - ocean_qflux) * dt / land_sea_heat_capacity
           delta_h_thermo_ice = 0
         elsewhere ! = ice-covered = ( h>0 )
           delta_t_ml = - ( thermo_ice_basal_flux_const * ( t_ml - t_thermo_ice_base ) - ocean_qflux ) * dt / land_sea_heat_capacity
           delta_h_thermo_ice = ( corrected_flux_noq - thermo_ice_basal_flux_const * ( t_ml - t_thermo_ice_base ) ) * dt / L_thermo_ice
         endwhere
         
         where ( t_ml + delta_t_ml .lt. t_surf_freeze ) ! = frazil growth =
           delta_h_thermo_ice = delta_h_thermo_ice - ( t_ml + delta_t_ml - t_surf_freeze ) * (land_sea_heat_capacity)/L_thermo_ice
           delta_t_ml = t_surf_freeze - t_ml
         endwhere
         
         where ( ( h_thermo_ice .gt. 0 ) .and. ( h_thermo_ice + delta_h_thermo_ice .le. 0 ) ) ! = complete ablation =
           delta_t_ml = delta_t_ml - ( h_thermo_ice + delta_h_thermo_ice ) * L_thermo_ice/land_sea_heat_capacity
           delta_h_thermo_ice = - h_thermo_ice
         endwhere
         
         ! = update surface temperature =
         ! compute surface melt
         where ( h_thermo_ice + delta_h_thermo_ice .gt. 0 ) ! surface is ice-covered
           ! calculate increment in steady-state ice surface temperature
           flux_thermo_ice = k_thermo_ice / (h_thermo_ice + delta_h_thermo_ice) * ( t_thermo_ice_base - t_surf )
           delta_t_thermo_ice = ( - corrected_flux_noq + flux_thermo_ice ) &
                         / ( k_thermo_ice / (h_thermo_ice + delta_h_thermo_ice) + dFdt_surf )
           where ( t_surf + delta_t_thermo_ice .gt. t_surf_freeze ) ! surface ablation
             delta_t_thermo_ice = t_surf_freeze - t_surf
           endwhere
           ! surface is ice-covered, so update t_surf as ice surface temperature
           delta_t_surf = delta_t_thermo_ice
           t_surf = t_surf + delta_t_thermo_ice
         elsewhere ! ice-free, so update t_surf as mixed layer temperature
         delta_t_surf = delta_t_ml
           t_surf = t_ml + delta_t_ml
         endwhere
         
   
   

        ! = update state =
        t_ml = t_ml + delta_t_ml
        h_thermo_ice = h_thermo_ice + delta_h_thermo_ice
        
       
    endwhere 

    if (coupled .and. ocean_tile) call coupled_thermo_ice()
    
    
    ! compute albedo
    where (land_ice_mask)
        albedo_out(:,:) = albedo_value * land_albedo_prefactor
    elsewhere 
        where ( h_thermo_ice .gt. 0.0 )
            albedo_out(:,:) = thermo_ice_albedo
        elsewhere 
            albedo_out(:,:) = albedo_value 
        endwhere
    endwhere 


    if (do_prescribe_albedo) then 
      if (prescribe_nh_albedo_only) then 
        where (rad_lat_2d .gt. 0) 
          albedo_out(:,:) = albedo_in(:,:)    
        endwhere 
      else 
        albedo_out(:,:) = albedo_in(:,:)
      endif 
    endif 
    
    if ( id_albedo > 0 .and. grid_diags) used = send_data ( id_albedo, albedo_out, Time_next ) 
endif 
  
if ((do_prescribe_albedo).and.(.not.do_thermo_ice)) then 

    albedo_out(:,:) = albedo_in(:,:)
    if ( id_albedo > 0 .and. grid_diags) used = send_data ( id_albedo, albedo_out, Time_next ) 
endif ! if (do_thermo_ice)

if (albedo_from_nudge_file) then 
  where ((h_thermo_ice .gt. 0) .and. (h_thermo_ice_in .le. 0)) 
    albedo_out(:,:) = albedo_value
  endwhere 
endif 

!
! Finally calculate the increments for the lowest atmospheric layer
!
if (do_explicit) then ! NTL 01/23 thermodynamic sea ice uses explicit step 
    Tri_surf%delta_t = fn_t 
    if (evaporation) Tri_surf%delta_tr(:,:,nhum) = fn_q 
else 
    Tri_surf%delta_t = fn_t + en_t * delta_t_surf
    if (evaporation) Tri_surf%delta_tr(:,:,nhum) = fn_q + en_q * delta_t_surf
endif 

if (coupled .and. surface_tiles) then
  ! the evaporation the atmosphere feels from this tile, for the water budget (see
  ! mixed_layer_tile_coupling)
  dts_flux = delta_t_surf
  if (do_explicit) dts_flux = 0.
  tile_evap(:,:,tile) = 0.
  if (evaporation) tile_evap(:,:,tile) = alpha_q + beta_q * dts_flux
endif

if (coupled .and. ice_tile) then
  ! fluxes over the ice, at the surface temperature the ice model sent, with their
  ! derivatives so it can apply them at each ice category's own temperature
  sens_out = CP_AIR * (alpha_t + beta_t * dts_flux)
  evap_out = 0.
  if (evaporation) evap_out = alpha_q + beta_q * dts_flux
  swnet_out = net_surf_sw_down   ! here the downward shortwave, before the ice albedo
  where (ice_alb_rad < 1.) swnet_out = net_surf_sw_down / (1. - ice_alb_rad)
  ! Precipitation falls as snow on the ice when the lowest-level air is below freezing
  ! (without snow SIS2's ice is bare, with a low albedo). Isca's condensation released
  ! only the latent heat of vaporisation, so the latent heat of fusion the snow carries
  ! is given back to the ice as heat, which keeps the total energy unchanged.
  snow_out = 0.
  if (present(t_atm)) then
    where (t_atm < TFREEZE_WATER) snow_out = precip
  endif
  sens_out = sens_out - HLF_WATER * snow_out
  call ocean_coupler_put_ice(Time, -flux_u, -flux_v, swnet_out, surf_lw_down, sens_out, &
                             dhdt_surf, evap_out, dedt_surf, precip - snow_out, snow_out)
  ice_alb_rad = albedo_out
endif

if (coupled .and. ocean_tile) then
  ! Send the fluxes consistent with the implicit surface step, i.e. the ones
  ! the atmosphere actually feels, so that energy and water are conserved.
  if (.not. (present(flux_u) .and. present(flux_v) .and. present(precip) .and. present(psurf))) &
    call error_mesg('mixed_layer', 'flux_u, flux_v, precip and psurf are needed for ocean coupling', FATAL)
  ! with explicit stepping (do_explicit, required by do_thermo_ice) the atmosphere feels the
  ! fluxes at the old surface temperature
  dts_flux = delta_t_surf
  if (do_explicit) dts_flux = 0.
  sens_out  = CP_AIR * (alpha_t + beta_t * dts_flux)
  lwnet_out = surf_lw_down - (alpha_lw + beta_lw * dts_flux)
  swnet_out = net_surf_sw_down
  if (evaporation) then
    evap_out = alpha_q + beta_q * dts_flux
  else
    evap_out = 0.
  endif
  ! Under coupled sea ice the ocean gets only the ice-ocean heat flux (basal flux and melt):
  ! the radiative and turbulent heat goes into the ice. Wind stress (free-drifting ice) and
  ! water (the zero-layer ice holds none) pass through; MOM6 takes HLV*evap as latent heat
  ! loss, so the sensible flux sent is adjusted to deliver exactly ice_ocean_heat.
  if (do_thermo_ice) then
    where (ice_covered_coupled)
      swnet_out = 0.
      lwnet_out = 0.
      sens_out = - ice_ocean_heat - HLV * evap_out
    end where
  endif
  ! flux_u/v from surface_flux are the stress on the atmosphere; the ocean feels minus that.
  ! Under prescribed sea ice (update_albedo_from_ice) the ocean gets no fluxes at all: the
  ! ice is treated as an insulating lid, as there is no sea-ice model in the coupled system.
  call ocean_coupler_put(Time, ice_lid(-flux_u), ice_lid(-flux_v), ice_lid(swnet_out), &
                         ice_lid(lwnet_out), ice_lid(sens_out), ice_lid(evap_out), &
                         ice_lid(precip), psurf)
  if (id_cpl_heat > 0) used = send_data(id_cpl_heat, &
        swnet_out + lwnet_out - sens_out - HLV*evap_out, Time)
  if (do_thermo_ice) then
    if (id_cpl_frazil > 0) used = send_data(id_cpl_frazil, frazil_coupled, Time)
    if (id_ice_ocean_heat > 0) used = send_data(id_ice_ocean_heat, ice_ocean_heat, Time)
  endif
  if (id_cpl_evap > 0) used = send_data(id_cpl_evap, evap_out, Time)
  if (id_cpl_precip > 0) used = send_data(id_cpl_precip, precip, Time)
  ! Water the ocean does not get directly becomes runoff. The precip/evap fields reach
  ! the ocean, through conservative remapping, from the ocean-covered fraction of each
  ! cell only, so the rest of the cell's P-E is routed (cells under the ice lid send
  ! nothing directly, so all of theirs). Land points are overridden by the bucket.
  uncoupled_pme = (precip - evap_out) * (1. - min(max(ocean_frac_coupled, 0.), 1.))
  if (update_albedo_from_ice) then
    where (ice_concentration > ice_concentration_threshold .and. .not. land_mask) &
      uncoupled_pme = precip - evap_out
  endif
endif

!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
! We have taken a time-step, send the values at the next time level.
if (grid_diags) then
  if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time_next)
  if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time_next)
  if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q_total, Time_next)
  if(id_delta_t_surf > 0)   used = send_data(id_delta_t_surf, delta_t_surf, Time_next)
endif
if (.not. ocean_tile) return   ! the rest describe the ocean (tile)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time_next)
! NTL 01/23 thermodynamic sea ice
if(id_h_thermo_ice > 0) used = send_data(id_h_thermo_ice, h_thermo_ice, Time_next)
if(id_t_ml > 0) used = send_data(id_t_ml, t_ml, Time_next)
if(id_flux_thermo_ice > 0) used = send_data(id_flux_thermo_ice, flux_thermo_ice, Time_next)

contains

subroutine coupled_thermo_ice()
! Thermodynamic sea ice over the points coupled to an external ocean. The zero-layer
! ice of the slab version above, with the ocean beneath it the coupled ocean model:
!  - open water: t_surf is the ocean's SST, all surface fluxes go to the ocean, and the
!    frazil heat the ocean had to add to stay at freezing grows ice instead;
!  - ice: the surface energy balance grows or melts the ice, which also takes the basal
!    flux thermo_ice_basal_flux_const*(SST - t_thermo_ice_base) from the ocean; ice that
!    melts completely returns the excess heat to the ocean.
! Heat per unit area of Isca cell for the ice, per unit area of ocean for the ocean: the
! ocean covers the fraction frac_cpl of a coupled cell, the ice the whole cell. Uses the
! explicit surface step (do_explicit), as the slab ice does.
h_ice_old = h_old_save
t_surf_old = t_surf_save
cpl = coupled_mask
frac_cpl = min(max(ocean_frac_coupled, 0.), 1.)
frazil_here = frazil_coupled
if (surface_tiles) then
  ! ice and ocean both cover the whole ocean tile; frazil arrives per unit cell area
  where (cpl) frazil_here = frazil_coupled / max(frac_cpl, 1.e-10)
  frac_cpl = 1.
endif
ice_covered_coupled = cpl .and. (h_ice_old > 0.)
basal_flux = 0. ; melt_excess = 0. ; ice_ocean_heat = 0.
where (ice_covered_coupled)
  basal_flux = thermo_ice_basal_flux_const * (sst_coupled - t_thermo_ice_base)
  delta_h_thermo_ice = (corrected_flux - basal_flux * frac_cpl) * dt / L_thermo_ice
elsewhere (cpl)
  delta_h_thermo_ice = 0.
end where
where (cpl) delta_h_thermo_ice = delta_h_thermo_ice + frazil_here * dt / L_thermo_ice
where (cpl .and. (h_ice_old + delta_h_thermo_ice < 0.))   ! complete ablation
  melt_excess = - (h_ice_old + delta_h_thermo_ice) * L_thermo_ice / dt
  delta_h_thermo_ice = - h_ice_old
end where
where (ice_covered_coupled) ice_ocean_heat = - basal_flux + melt_excess / max(frac_cpl, 1.e-3)
where (cpl)
  h_thermo_ice = h_ice_old + delta_h_thermo_ice
  t_ml = sst_coupled
end where
! surface temperature: steady-state ice surface where ice remains, else the SST
where (ice_covered_coupled .and. h_thermo_ice > 0.)
  flux_thermo_ice = k_thermo_ice / h_thermo_ice * (t_thermo_ice_base - t_surf_old)
  delta_t_thermo_ice = (- corrected_flux + flux_thermo_ice) / (k_thermo_ice / h_thermo_ice + dFdt_surf)
  where (t_surf_old + delta_t_thermo_ice > t_surf_freeze) delta_t_thermo_ice = t_surf_freeze - t_surf_old
  t_surf = t_surf_old + delta_t_thermo_ice
elsewhere (cpl)
  flux_thermo_ice = 0.
  t_surf = sst_coupled
end where
where (cpl) delta_t_surf = t_surf - t_surf_old
! Ocean can lie under a cell that Isca treats as uncoupled (a polar cap, or land where the
! coastlines differ). Its frazil grows that cell's slab ice, so the heat is not lost.
! Under Isca land it is still dropped until Isca has fractional surface tiles.
where (.not. cpl .and. .not. land_mask) &
  h_thermo_ice = h_thermo_ice + frazil_coupled * dt / L_thermo_ice
end subroutine coupled_thermo_ice

end subroutine mixed_layer

!=================================================================================================================================

function ice_lid(field)
! field, set to zero over ocean points covered by prescribed sea ice
real, intent(in), dimension(:,:) :: field
real, dimension(size(field,1),size(field,2)) :: ice_lid
ice_lid = field
if (update_albedo_from_ice) then
  where (ice_concentration > ice_concentration_threshold .and. .not. land_mask) ice_lid = 0.
endif
end function ice_lid

!=================================================================================================================================

subroutine mixed_layer_uncoupled_pme(pme)
! P-E [kg/m2/s] over every point not covered by the coupled ocean, from the
! last mixed_layer call (zero when not coupled). The caller routes it to the
! ocean, except where a land bucket stores the water instead.
real, intent(out), dimension(:,:) :: pme
if (allocated(uncoupled_pme)) then
  pme = uncoupled_pme
else
  pme = 0.
endif
end subroutine mixed_layer_uncoupled_pme

!=================================================================================================================================

subroutine albedo_calc(albedo_inout,Time)

real, intent(out), dimension(:,:) :: albedo_inout
type(time_type), intent(in)       :: Time

albedo_inout=albedo_initial

if(update_albedo_from_ice) then

  if(trim(ice_albedo_method) == 'step_function') then
    where(ice_concentration.gt.ice_concentration_threshold)
      albedo_inout=ice_albedo_value
    end where
  else if(trim(ice_albedo_method) == 'ramp_function') then
    albedo_inout = albedo_inout*(1.0-ice_concentration) + ice_albedo_value*ice_concentration
  else
    call error_mesg('mixed_layer','"'//trim(ice_albedo_method)//'"'//' is not a valid method for determining'// &
      'albedo when ice is present. Choices are: step_function or ramp_function.', FATAL)
  endif

    if ( id_ice_conc > 0 ) used = send_data ( id_ice_conc, ice_concentration, Time )
    if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo_inout, Time )

endif

end subroutine albedo_calc
!=================================================================================================================================

subroutine read_ice_conc(Time)

type(time_type), intent(in)       :: Time


call interpolator( ice_interp, Time, ice_concentration, trim(ice_file_name) )
if ( id_ice_conc > 0 ) used = send_data ( id_ice_conc, ice_concentration, Time )

end subroutine read_ice_conc
!=================================================================================================================================

subroutine mixed_layer_tile_albedo(albedo, albedo_tile)
! Initial albedo of each surface tile (surface_tiles). With do_thermo_ice or
! do_prescribe_albedo the albedo is part of the state, and all tiles start from the
! albedo of the cell.
real, intent(in),  dimension(:,:)   :: albedo
real, intent(out), dimension(:,:,:) :: albedo_tile
integer :: k
if (do_thermo_ice .or. do_prescribe_albedo) then
  do k = 1, size(albedo_tile, 3)
    albedo_tile(:,:,k) = albedo
  enddo
else
  albedo_tile = albedo_initial_tile
endif
end subroutine mixed_layer_tile_albedo

!=================================================================================================================================

subroutine mixed_layer_tile_coupling(ocean_frac, evap, ice_frac)
! With surface_tiles and a coupled ocean: the fraction of each cell the ocean model covers,
! as last received, each tile's evaporation from the last step [kg/m2/s], and the fraction
! of that ocean covered by the ocean component's sea ice (zero without a sea-ice model).
real, intent(out), dimension(:,:)   :: ocean_frac, ice_frac
real, intent(out), dimension(:,:,:) :: evap
ocean_frac = 0.
evap = 0.
ice_frac = 0.
if (.not. coupled) return
ocean_frac = min(max(ocean_frac_coupled, 0.), 1.)
if (allocated(tile_evap)) evap = tile_evap
if (coupled_ice) ice_frac = min(max(ice_frac_cpl, 0.), 1.)
end subroutine mixed_layer_tile_coupling

!=================================================================================================================================

subroutine mixed_layer_tile_diagnostics(Time_next, t_surf, flux_t, flux_lhe, delta_t_surf, albedo)
! With surface_tiles: the mixed_layer diagnostics that describe the whole cell, as
! averages over the tiles weighted by area.
type(time_type), intent(in) :: Time_next
real, intent(in), dimension(:,:) :: t_surf, flux_t, flux_lhe, delta_t_surf, albedo
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time_next)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time_next)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, flux_lhe, Time_next)
if(id_delta_t_surf > 0) used = send_data(id_delta_t_surf, delta_t_surf, Time_next)
if(id_albedo > 0 .and. (do_thermo_ice .or. do_prescribe_albedo)) used = send_data(id_albedo, albedo, Time_next)
end subroutine mixed_layer_tile_diagnostics

!=================================================================================================================================

subroutine mixed_layer_end(t_surf, &
                           h_thermo_ice, t_ml, const_correct, nudge_out, albedo, & ! NTL 01/23 thermodynamic sea ice
                           bucket_depth, restart_file_bucket_depth) 

real, intent(inout), dimension(:,:) :: t_surf
real, intent(inout), dimension(:,:) :: h_thermo_ice, t_ml, const_correct, nudge_out, albedo ! NTL 01/23 thermodynamic sea ice
real, intent(inout), dimension(:,:,:) :: bucket_depth
logical, intent(in)                 :: restart_file_bucket_depth
integer:: unit

if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
if (do_thermo_ice) then ! NTL 01/23 thermodynamic sea ice
    call write_data(trim('RESTART/mixed_layer.res'), 'h_thermo_ice',   h_thermo_ice, grid_domain)
    call write_data(trim('RESTART/mixed_layer.res'), 't_ml',   t_ml, grid_domain)
    call write_data(trim('RESTART/mixed_layer.res'), 'const_correct',   const_correct, grid_domain)
    call write_data(trim('RESTART/mixed_layer.res'), 'nudge_out',   nudge_out, grid_domain)
    call write_data(trim('RESTART/mixed_layer.res'), 'albedo',   albedo, grid_domain)
else if (do_prescribe_albedo) then 
    call write_data(trim('RESTART/mixed_layer.res'), 'albedo',   albedo, grid_domain)
endif 
if (restart_file_bucket_depth) then
    call write_data(trim('RESTART/mixed_layer.res'), 'bucket_depth',   bucket_depth, grid_domain)
endif

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
