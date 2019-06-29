                  module physics_driver_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <OVERVIEW>
!     Provides high level interfaces for calling the entire
!     FMS atmospheric physics package.
!
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
! </OVERVIEW>
! <DESCRIPTION>
!     This version of physics_driver_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Radiation, Rayleigh damping, gravity wave drag, vertical diffusion of
!     momentum and tracers, and the downward pass of vertical diffusion for
!     temperature and specific humidity are performed in the down routine.
!     The up routine finishes the vertical diffusion and computes moisture
!     related terms (convection,large-scale condensation, and precipitation).
! </DESCRIPTION>
! <DIAGFIELDS>
! </DIAGFIELDS>
! <DATASET NAME="physics_driver.res">
! native format restart file
! </DATASET>
!
! <DATASET NAME="physics_driver.res.nc">
! netcdf format restart file
! </DATASET>


! <INFO>

!   <REFERENCE>            </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!   </NOTE>
!   <FUTURE> Deal with conservation of total energy?              </FUTURE>

! </INFO>
!   shared modules:

use time_manager_mod,        only: time_type, get_time, operator (-), &
                                   time_manager_init
use field_manager_mod,       only: field_manager_init, MODEL_ATMOS
use tracer_manager_mod,      only: tracer_manager_init, &
                                   get_number_tracers, &
                                   get_tracer_names

! use atmos_tracer_driver_mod, only: atmos_tracer_driver_init,    &
!                                    atmos_tracer_driver_time_vary, &
!                                    atmos_tracer_driver_endts, &
!                                    atmos_tracer_driver,  &
!                                    atmos_tracer_driver_end
use mpp_mod,                 only: input_nml_file
use fms_mod,                 only: mpp_clock_id, mpp_clock_begin,   &
                                   mpp_clock_end, CLOCK_MODULE_DRIVER, &
                                   fms_init,  &
                                   open_namelist_file, stdlog, stdout,  &
                                   write_version_number, field_size, &
                                   file_exist, error_mesg, FATAL,   &
                                   WARNING, NOTE, check_nml_error, &
                                   close_file, mpp_pe, mpp_root_pe, &
                                   mpp_error, mpp_chksum
use fms_io_mod,              only: restore_state, &
                                   register_restart_field, restart_file_type, &
                                   save_restart, get_mosaic_tile_file
use constants_mod,           only: RDGAS, PI

use diag_manager_mod,        only: register_diag_field, send_data

!    shared radiation package modules:

use rad_utilities_mod,       only: aerosol_type, radiative_gases_type, &
                                   rad_utilities_init, rad_output_type,&
                                   cld_specification_type,   &
                                   surface_type, &
                                   atmos_input_type, microphysics_type

!    component modules:

use  moist_processes_mod,    only: moist_processes,    &
                                   moist_processes_init,  &
                                   set_cosp_precip_sources, &
                                   moist_processes_endts, &
                                   moist_processes_end,  &
                                   moist_alloc_init, &
                                   moist_alloc_end, &
                                   doing_strat!,          &
!                                   moist_processes_restart
! use moistproc_kernels_mod,   only:  moistproc_init, moistproc_end

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end, &
                                   vert_turb_driver_restart

use idealized_moist_phys_mod, only: idealized_moist_phys_init , idealized_moist_phys , idealized_convection_and_lscale_cond, idealized_radiation_and_optional_surface_flux, idealized_moist_phys_end, surf_diff_type

! use atmosphere_mod,           only: get_nhum

use         mpp_domains_mod, only: domain2D !s added to enable land reading

use vert_diff_driver_mod,    only: vert_diff_driver_down,  &
                                   vert_diff_driver_up, &
                                   vert_diff_driver_init
                                !    vert_diff_driver_end,   &
                                !    surf_diff_type

! use radiation_driver_mod,    only: radiation_driver_init,    &
!                                    define_rad_times, define_surface,   &
!                                    define_atmos_input_fields, &
!                                    radiation_driver_time_vary, &
!                                    radiation_driver_endts, &
!                                    radiation_driver,  &
!                                    return_cosp_inputs, &
!                                    atmos_input_dealloc,    &
!                                    microphys_dealloc, &
!                                    surface_dealloc, &
!                                    radiation_driver_end, &
!                                    radiation_driver_restart
  
! use cloud_spec_mod,          only: cloud_spec_init, cloud_spec, &
!                                    cloud_spec_dealloc, cloud_spec_end

! use aerosol_mod,             only: aerosol_init, aerosol_driver, &
!                                    aerosol_time_vary, &
!                                    aerosol_endts, &
!                                    aerosol_dealloc, aerosol_end
 
! use radiative_gases_mod,     only: radiative_gases_init,   &
!                                    radiative_gases_time_vary, &
!                                    radiative_gases_endts, &
!                                    define_radiative_gases, &
!                                    radiative_gases_dealloc, &
!                                    radiative_gases_end,     &
!                                    radiative_gases_restart

use damping_driver_mod,      only: damping_driver, &
!                                    damping_driver_init, &
!                                    damping_driver_endts, &
                                   damping_driver_end
!                                    damping_driver_restart

! use grey_radiation_mod,       only: grey_radiation_init, grey_radiation, &
!                                     grey_radiation_end

!--> h1g, cjg
! use clubb_driver_mod,         only: clubb_init, clubb, clubb_end
!<-- h1g, cjg

#ifdef SCM
! Option to add SCM radiative tendencies from forcing to lw_tendency
! and radturbten

use scm_forc_mod,            only: use_scm_rad, add_scm_tdtlw, add_scm_tdtsw

#endif

!-----------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id: physics_driver.F90,v 20.0 2013/12/13 23:18:40 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public  physics_driver_init, physics_driver_down,   &
        ! physics_driver_down_endts, physics_driver_up_endts, &
        physics_driver_up, physics_driver_end, &
        do_moist_in_phys_up, get_diff_t, &
        get_radturbten, zero_radturbten, physics_driver_restart

private          &

!  called from physics_driver_down:
         check_args, &

!  called from check_args:
         check_dim

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface


!---------------------------------------------------------------------
!------- namelist ------

logical :: do_moist_processes = .true.
                               ! call moist_processes routines
real    :: tau_diff = 3600.    ! time scale for smoothing diffusion 
                               ! coefficients
integer :: do_clubb = 0        ! activate clubb parameterization ?
logical :: do_cosp = .false.   ! activate COSP simulator ?
logical :: do_modis_yim = .true. ! activate simple modis simulator ?
logical :: do_radiation = .true.
                               ! calculating radiative fluxes and
                               ! heating rates?
logical :: do_grey_radiation = .false. ! do grey radiation scheme?
real    :: R1 = 0.25           ! rif:(09/10/09) In Grey radiation we are computing just the total   
real    :: R2 = 0.25           ! SW radiation. We need to divide it into 4 components
real    :: R3 = 0.25           ! to go through the Coupler and Ice modules.
real    :: R4 = 0.25           ! 	Sum[R(i)*SW] = SW  


real    :: diff_min = 1.e-3    ! minimum value of a diffusion 
                               ! coefficient beneath which the
                               ! coefficient is reset to zero
logical :: diffusion_smooth = .true.
                               ! diffusion coefficients should be 
                               ! smoothed in time?
logical :: use_cloud_tracers_in_radiation = .true.
                               ! if true, use lsc cloud tracer fields
                               ! in radiation (these transported on
                               ! current step, will have non-realizable
                               ! total cloud areas at some points); if
                               ! false, then use balanced (realizable)
                               ! fields saved at end of last step
                               ! only an issue when both lsc and conv
                               ! clouds are active (AM3)
logical :: donner_meso_is_largescale = .true.
                               ! donner meso clouds are treated as 
                               ! largescale (rather than convective)
                               ! as far as the COSP simulator is 
                               ! concerned ?
logical :: allow_cosp_precip_wo_clouds = .true.
                               ! COSP will see {ls, cv} precip in grid-
                               ! boxes w/o {ls, cv} clouds ?
logical :: override_aerosols_radiation = .false.
                               ! use offline aerosols for radiation 
                               ! calculation
                               ! (via data_override in aerosol_driver)?
logical :: override_aerosols_cloud = .false.
                               ! use offline aerosols for cloud calculation
                               ! (via data_override in aerosol_driver)?
character(len=16) :: cosp_precip_sources = '    '
                               ! sources of the precip fields to be sent to
                               ! COSP. Default = '  ' implies precip from
                               ! the radiatively-active clouds defined by
                               ! variable cloud_type_form in cloud_spec_nml
                               ! will be sent. Other available choices:
                               ! 'strat', 'deep', 'uw', 'stratdeep',
                               ! 'stratuw', deepuw', 'stratdeepuw', 
                               ! 'noprecip'. 
                               ! CURRENTLY NOT AVAILABLE: precip from ras,
                               ! lsc and mca. For completeness, these could
                               ! be made available, but since no cloud 
                               ! fields are saved for these schemes to be
                               ! made available to COSP, they are also
                               ! currently not considered as precip
                               ! sources.

! <NAMELIST NAME="physics_driver_nml">
!  <DATA NAME="do_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!calculating radiative fluxes and
! heating rates?
!  </DATA>
!  <DATA NAME="do_moist_processes" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!call moist_processes routines
!  </DATA>
!  <DATA NAME="tau_diff" UNITS="" TYPE="real" DIM="" DEFAULT="3600.">
!time scale for smoothing diffusion 
! coefficients
!  </DATA>
!  <DATA NAME="diff_min" UNITS="" TYPE="real" DIM="" DEFAULT="1.e-3">
!minimum value of a diffusion 
! coefficient beneath which the
! coefficient is reset to zero
!  </DATA>
!  <DATA NAME="diffusion_smooth" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!diffusion coefficients should be 
! smoothed in time?
!  </DATA>
! <DATA NAME="override_aerosols_radiation" UNITS="" TYPE="logical" DIM=""  DEFAULT=".false.">
!use offline aerosols for radiation calculation
! (via data_override in aerosol_driver)?
!  </DATA>
!  <DATA NAME="override_aerosols_cloud" UNITS="" TYPE="logical" DIM="" DEFAULT=".false.">
!use offline aerosols for cloud calculation
! (via data_override in aerosol_driver)?
!  </DATA>
! </NAMELIST>
!

! ---> h1g, 2012-08-28, add option of applying surface fluxes in host-model
!                       by default .true. (that is, applying surface fluxes in host-model)
logical :: l_host_applies_sfc_fluxes = .true.
logical :: print_s_messages = .true.
! <--- h1g, 2012-08-28

namelist / physics_driver_nml / do_radiation, &
                                do_clubb,  &  ! cjg, h1g
                                do_cosp, &
                                do_modis_yim, &
                                donner_meso_is_largescale, &
                                allow_cosp_precip_wo_clouds, &
                                do_moist_processes, tau_diff,      &
                                diff_min, diffusion_smooth, &
                                use_cloud_tracers_in_radiation, &
                                R1, R2, R3, R4,  &
                                override_aerosols_radiation,  &
                                override_aerosols_cloud,    &
                                cosp_precip_sources,    &
                                l_host_applies_sfc_fluxes, &
                                print_s_messages

!---------------------------------------------------------------------
!------- public data ------
! <DATA NAME="surf_diff_type" UNITS="" TYPE="surf_diff_type" DIM="" DEFAULT="">
! Defined in vert_diff_driver_mod, republished here. See vert_diff_mod for details.
! </DATA>

public  surf_diff_type   ! defined in  vert_diff_driver_mod, republished
                         ! here
 
!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
! list of restart versions readable by this module:
!
! version 1: initial implementation 1/2003, contains diffusion coef-
!            ficient contribution from cu_mo_trans_mod. This variable
!            is generated in physics_driver_up (moist_processes) and
!            used on the next step in vert_diff_down, necessitating
!            its storage.
!
! version 2: adds pbltop as generated in vert_turb_driver_mod. This 
!            variable is then used on the next timestep by topo_drag
!            (called from damping_driver_mod), necessitating its 
!            storage.
!
! version 3: adds the diffusion coefficients which are passed to 
!            vert_diff_driver.  These diffusion are saved should
!            smoothing of vertical diffusion coefficients be turned
!            on.
!
! version 4: adds a logical variable, convect, which indicates whether
!            or not the grid column is convecting. This diagnostic is
!            needed by the entrain_module in vert_turb_driver.
!
! version 5: adds radturbten when strat_cloud_mod is active, adds 
!            lw_tendency when edt_mod or entrain_mod is active.
!
! version 6: adds donner cell and meso cloud variables when donner_deep
!            is activated.

! version 7: adds shallow convection cloud variables when uw_conv
!            is activated.

! version 8: adds lsc cloud props for radiation. only readable when in
!            netcdf mode.


!---------------------------------------------------------------------
integer, dimension(8) :: restart_versions = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)

!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics 
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!  
!    diff_cu_mo     contains contribution to difusion coefficient
!                   coming from cu_mo_trans_mod (called from 
!                   moist_processes in physics_driver_up) and then used 
!                   as input on the next time step to vert_diff_down 
!                   called in physics_driver_down.
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    radturbten     the sum of the radiational and turbulent heating,
!                   generated in both physics_driver_down (radiation)
!                   and physics_driver_up (turbulence) and then used
!                   in moist_processes
!    lw_tendency    longwave heating rate, generated in radiation and
!                   needed in vert_turb_driver when either edt_mod
!                   or entrain_mod is active. must be saved because
!                   radiation is not calculated on each step.
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down        
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from 
!                   physics_driver_down on the next step.
!----------------------------------------------------------------------
real,    dimension(:,:,:), allocatable :: diff_cu_mo, diff_t, diff_m
real,    dimension(:,:,:), allocatable :: radturbten, lw_tendency
real,    dimension(:,:)  , allocatable :: pbltop, cush, cbmf
logical, dimension(:,:)  , allocatable :: convect
real,    dimension(:,:,:), allocatable ::       &
                           cell_cld_frac, cell_liq_amt, &
                           cell_liq_size, cell_ice_amt, cell_ice_size, &
                           cell_droplet_number, &
                           meso_cld_frac, meso_liq_amt, meso_liq_size, &
                           meso_ice_amt, meso_ice_size, &
                           meso_droplet_number, &
                           lsc_cloud_area, lsc_liquid, &
                           lsc_ice, lsc_droplet_number, &
                           lsc_ice_number, &
             !snow, rain in radiation
                           lsc_rain, lsc_snow,              &
                           lsc_rain_size, lsc_snow_size,              &
                           shallow_cloud_area, shallow_liquid, &
                           shallow_ice, shallow_droplet_number, &
!cms++ not yet activated, but may be needed ...
                           shallow_ice_number, &
!cms--
                           temp_last, q_last
real,    dimension(:,:,:,:), allocatable ::  tau_stoch, lwem_stoch, &
                           stoch_cloud_type, stoch_conc_drop, &
                           stoch_conc_ice, stoch_size_drop, &
                           stoch_size_ice
real,    dimension(:,:,:), allocatable ::  fl_lsrain, fl_lssnow, &
                                           fl_lsgrpl, &
                                           fl_donmca_rain, fl_donmca_snow,&
                                           fl_ccrain, fl_ccsnow, &
                                           mr_ozone
real,       dimension(:,:), allocatable  :: daytime
integer,    dimension(:,:)  , allocatable :: nsum_out
real   ,    dimension(:,:)  , allocatable :: tsurf_save

! --->h1g
real,    dimension(:,:,:), allocatable ::  dcond_ls_liquid, dcond_ls_ice
real,    dimension(:,:,:), allocatable ::  Ndrop_act_CLUBB,  Icedrop_act_CLUBB
real,    dimension(:,:,:), allocatable ::  ndust, rbar_dust
real,    dimension(:,:,:), allocatable ::  diff_t_clubb
! <---h1g
   
!--- for netcdf restart
type(restart_file_type), pointer, save :: Phy_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer                                :: vers
integer                                :: now_doing_strat  
integer                                :: now_doing_entrain
integer                                :: now_doing_edt
real, allocatable                      :: r_convect(:,:)

!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
integer :: radiation_clock, damping_clock, turb_clock,   &
           tracer_clock, diff_up_clock, diff_down_clock, &
           moist_processes_clock, cosp_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
logical   :: do_check_args = .true.   ! argument dimensions should 
                                      ! be checked ?
logical   :: module_is_initialized = .false.
                                      ! module has been initialized ?
logical   :: doing_edt                ! edt_mod has been activated ?
logical   :: doing_entrain            ! entrain_mod has been activated ?
logical   :: doing_donner             ! donner_deep_mod has been 
                                      ! activated ?
logical   :: doing_uw_conv            ! uw_conv shallow cu mod has been 
                                      ! activated ?
logical   :: doing_liq_num = .false.  ! Prognostic cloud droplet number has 
                                      ! been activated?
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
integer   :: ncol                     ! number of stochastic columns
 
type(radiative_gases_type), save   :: Rad_gases_tv
type(time_type)  :: Rad_time
logical          ::    need_aerosols, need_clouds, need_gases,   &
                       need_basic
logical    :: do_strat
integer   :: num_uw_tracers


!---------------------------------------------------------------------
!---------------------------------------------------------------------

character(len=4)     :: mod_name = 'phys'
character(len=32)    :: tracer_units, tracer_name
  character(len=128) :: diaglname
real                 :: missing_value = -999.
logical              :: step_to_call_cosp = .false.
logical              :: include_donmca_in_cosp

!integer                            :: id_tdt_phys, id_qdt_phys, &   ! cjg
integer                            :: id_tdt_phys,         &
                                      id_tdt_phys_vdif_dn, &
                                      id_tdt_phys_vdif_up, &
                                      id_tdt_phys_turb,    &
                                      id_tdt_phys_moist

integer, dimension(:), allocatable :: id_tracer_phys,         &  ! cjg
                                      id_tracer_phys_vdif_dn, &
                                      id_tracer_phys_vdif_up, &
                                      id_tracer_phys_turb,    &
                                      id_tracer_phys_moist


                            contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="physics_driver_init">
!  <OVERVIEW>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_init (Time, lonb, latb, axes, pref, &
!                             trs, Surf_diff, phalf, mask, kbot  )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!  </IN>
!  <INOUT NAME="trs" TYPE="real">
!   atmospheric tracer fields
!  </INOUT>
!  <INOUT NAME="Surf_diff" TYPE="surf_diff_type">
!   surface diffusion derived type
!  </INOUT>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model interface levels
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! <ERROR MSG="physics_driver_init must be called first" STATUS="FATAL">
! </ERROR>
! </SUBROUTINE>
!
!-->cjg
!subroutine physics_driver_init (Time, lonb, latb, axes, pref, &
subroutine physics_driver_init (Time, lonb, latb, lon, lat, axes, pref, &
                                trs, Surf_diff, phalf, grid_domain_in, mask, kbot, &
                                diffm, difft, is_in, ie_in, js_in, je_in, num_levels_in, nhum_in, &
                                surf_geopotential, Time_step  )
!<--cjg
!---------------------------------------------------------------------
!    physics_driver_init is the constructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type),         intent(in)              :: Time
real,dimension(:,:),     intent(in)              :: lonb, latb
real,dimension(:,:),     intent(in)              :: lon, lat    ! cjg
integer,dimension(4),    intent(in)              :: axes
real,dimension(:,:),     intent(in)              :: pref
real,dimension(:,:,:,:), intent(inout)           :: trs
type(surf_diff_type),    intent(inout)           :: Surf_diff
real,dimension(:,:,:),   intent(in)              :: phalf
type(domain2D),          intent(in)              :: grid_domain_in
real,dimension(:,:,:),   intent(in),   optional  :: mask
integer,dimension(:,:),  intent(in),   optional  :: kbot
real, dimension(:,:,:),  intent(out),  optional  :: diffm, difft

integer,                 intent(in),   optional  :: is_in, ie_in, js_in, je_in, num_levels_in, nhum_in
real, dimension(:,:),    intent(in),   optional  :: surf_geopotential
type(time_type),         intent(in),   optional  :: Time_step

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     Time       current time (time_type)
!     lonb       longitude of the grid box corners [ radians ]
!     latb       latitude of the grid box corners [ radians ]
!     axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!     pref       two reference profiles of pressure at nlev+1 levels
!                pref(nlev+1,1)=101325. and pref(nlev+1,2)=81060.
!     phalf      pressure at model interface levels
!                [ Pa ]
!
!   intent(inout) variables:
!
!     trs        atmosperic tracer fields
!     Surf_diff  surface diffusion derived type variable
!
!   intent(in), optional variables:
!
!        mask    present when running eta vertical coordinate,
!                mask to remove points below ground
!        kbot    present when running eta vertical coordinate,
!                index of lowest model level above ground
!   
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(lonb,1)-1, size(latb,2)-1) :: sgsmtn
      character(len=64), dimension(:), pointer :: aerosol_names => NULL()
      character(len=64), dimension(:), pointer :: aerosol_family_names => NULL()
      integer          ::  id, jd, kd, n
      integer          ::  ierr, io, unit, logunit, outunit
      integer          ::  ndum
      character(len=16)::  cloud_type_form_out  ! indicator of radiatively 
                                                ! active clouds
      character(len=16)::  cosp_precip_sources_modified

      integer          ::  moist_processes_init_clock, damping_init_clock, &
                           turb_init_clock, diff_init_clock, &
                           cloud_spec_init_clock, aerosol_init_clock, &
                           grey_radiation_init_clock , radiative_gases_init_clock, &
                           radiation_init_clock, tracer_init_clock, &
                           cosp_init_clock

      real, dimension(size(lonb,1), size(lonb,2)) :: rad_lonb, rad_latb
      real, dimension(size(lon,1), size(lon,2))   :: rad_lon, rad_lat

!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       aerosol_names names associated with the activated aerosols
!                     that will be seen by the radiation package
!       aerosol_family_names
!              names associated with the activated aerosol
!              families that will be seen by the radiation package
!       id,jd,kd      model dimensions on the processor  
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called 
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init
      call tracer_manager_init
      call field_manager_init (ndum)
 
!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=physics_driver_nml, iostat=io)
      ierr = check_nml_error(io,"physics_driver_nml")
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=physics_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'physics_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!       if(do_radiation .and. do_grey_radiation) & 
!         call error_mesg('physics_driver_init','do_radiation and do_grey_radiation cannot both be .true.',FATAL)

!CONVERT INCOMING DIMENSIONAL ARRAYS INTO RADIANS

rad_lat = lat 
rad_latb = latb 
rad_lon = lon 
rad_lonb = lonb 

  call idealized_moist_phys_init(is_in, ie_in, js_in, je_in, num_levels_in, axes, surf_geopotential, Time, Time_step, nhum_in, rad_lon, rad_lat, rad_lonb, rad_latb, grid_domain_in, Surf_diff, do_grey_radiation=do_grey_radiation)

  id = size(lonb,1)-1 
  jd = size(latb,2)-1 
  kd = num_levels_in
  call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
  num_prog=ntp)

  call  moist_processes_init (id, jd, kd, lonb, latb, lon, lat, phalf, pref(:,1),&
  axes, Time, doing_donner)

  if (print_s_messages) write(6,*) 'checking nt ntp', nt, ntp
  call moist_alloc_init (id,jd,kd,nt,ntp) ! Are these the right way around? nt and ntp correct?

  call vert_diff_driver_init (Surf_diff, id, jd, kd, axes, Time, do_clubb )



!---------------------------------------------------------------------
!    allocate space for the module variables.
!---------------------------------------------------------------------
  allocate ( diff_t     (id, jd, kd) ) ; diff_t = 0.0
  allocate ( diff_m     (id, jd, kd) ) ; diff_m = 0.0
  allocate ( diff_cu_mo (id, jd, kd) ) ; diff_cu_mo = 0.0
  allocate ( pbltop     (id, jd) )     ; pbltop     = -999.0
  allocate ( cush       (id, jd) )     ; cush=-1. !miz
  allocate ( cbmf       (id, jd) )     ; cbmf=0.0 !miz
  allocate ( convect    (id, jd) )     ; convect = .false.
  allocate ( radturbten (id, jd, kd))  ; radturbten = 0.0
  allocate ( lw_tendency(id, jd, kd))  ; lw_tendency = 0.0
  allocate ( r_convect  (id, jd) )     ; r_convect   = 0.0


  allocate (fl_lsrain  (id, jd, kd))
  allocate (fl_lssnow  (id, jd, kd))
  allocate (fl_lsgrpl  (id, jd, kd))
  allocate (fl_ccrain  (id, jd, kd))
  allocate (fl_ccsnow  (id, jd, kd))
  allocate (fl_donmca_snow  (id, jd, kd))
  allocate (fl_donmca_rain  (id, jd, kd))
  allocate (mr_ozone   (id, jd, kd))
  allocate (daytime    (id, jd    ))
  allocate ( temp_last (id, jd, kd))
  allocate ( q_last    (id, jd, kd))
  fl_lsrain = 0.
  fl_lssnow = 0.
  fl_lsgrpl = 0.
  fl_ccrain = 0.
  fl_ccsnow = 0.
  fl_donmca_rain = 0.
  fl_donmca_snow = 0.
  mr_ozone  = 0.
  daytime =   0.
  temp_last = 0.
  q_last    = 0. 



!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
               write(logunit, nml=physics_driver_nml)
 
!---------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------


!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------



 end subroutine physics_driver_init


!######################################################################
! <SUBROUTINE NAME="physics_driver_down_time_vary">
!  <OVERVIEW>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down_time_vary (Time, Time_next)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time of next time step
!  </IN>
! </SUBROUTINE>
!

! subroutine physics_driver_down_time_vary (Time, Time_next, gavg_rrv, dt)

! !---------------------------------------------------------------------
! !    physics_driver_down_time_vary makes sure that all time-dependent, 
! !    spacially-independent calculations are completed before entering window 
! !    or thread loops. Resultant fields are usually saved as module variables in 
! !    the module where needed.
! !-----------------------------------------------------------------------

! type(time_type),         intent(in)             :: Time, Time_next
! real, dimension(:),      intent(in)             :: gavg_rrv
! real,                    intent(in)             :: dt

! !---------------------------------------------------------------------      
!       if (do_radiation) then
! !----------------------------------------------------------------------
! !    call define_rad_times to obtain the time to be used in the rad-
! !    iation calculation (Rad_time) and to determine which, if any, 
! !    externally-supplied inputs to radiation_driver must be obtained on 
! !    this timestep.  logical flags are returned indicating the need or 
! !    lack of need for the aerosol fields, the cloud fields, the rad-
! !    iative gas fields, and the basic atmospheric variable fields.
! !----------------------------------------------------------------------
!         call define_rad_times (Time, Time_next, Rad_time, &
!                                need_aerosols, need_clouds, &
!                                need_gases, need_basic)
!         call aerosol_time_vary (Rad_time)
!         call radiative_gases_time_vary (Rad_time, gavg_rrv,  &
!                                                            Rad_gases_tv)
!         call radiation_driver_time_vary (Rad_time, Rad_gases_tv)
        
! !--------------------------------------------------------------------
! !    define step_to_call_cosp to indicate that this is a radiation
! !    step and therefore one on which COSP should be called in 
! !    physics_driver_up.
! !--------------------------------------------------------------------
!         if (need_basic) then
!           step_to_call_cosp = .true.
!         else
!           step_to_call_cosp = .false.
!         endif

!       endif
!       call damping_driver_time_vary (dt)
!       call atmos_tracer_driver_time_vary (Time)


! !-------------------------------------------------------------------------      

! end subroutine physics_driver_down_time_vary



!######################################################################

! subroutine physics_driver_down_endts(is,js)

! integer, intent(in)  :: is,js

!       call damping_driver_endts
!       call atmos_tracer_driver_endts

!       IF (do_radiation) THEN
!         !  CALL aerosol_endts
!          CALL radiation_driver_endts  (is, js, Rad_gases_tv)
!          CALL radiative_gases_endts
!       END IF
   

! !--------------------------------------------------------------------
! !    set a flag to indicate that this check was done and need not be
! !    done again.
! !--------------------------------------------------------------------
!       do_check_args = .false.


! end subroutine physics_driver_down_endts

!######################################################################

subroutine physics_driver_up_endts (is,js)

integer, intent(in)  :: is,js

      call moist_processes_endts (is,js)
      ! call aerosol_endts

end subroutine physics_driver_up_endts


! !#####################################################################

! !--> cjg: code modification to allow diagnostic tracers in physics_up (20120508)
! !         
! !subroutine physics_driver_moist_init (ix,jx,kx,lx)
! !
! !integer, intent(in) :: ix,jx, kx, lx 
! !
! !      call moist_alloc_init (ix,jx,kx,lx)

! subroutine physics_driver_moist_init (ix,jx,kx,lx,mx)


! integer, intent(in) :: ix,jx, kx, lx, mx


!       call moist_alloc_init (ix,jx,kx,lx,mx)
! !<--cjg
!       call moistproc_init (ix,jx,kx, num_uw_tracers, do_strat)

! end subroutine physics_driver_moist_init 



! !######################################################################

! subroutine physics_driver_moist_end                  

!       call moist_alloc_end
!       call moistproc_end (do_strat)


! end subroutine physics_driver_moist_end                  




!######################################################################
! <SUBROUTINE NAME="physics_driver_down">
!  <OVERVIEW>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.    
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down (is, ie, js, je,                       &
!                                Time_prev, Time, Time_next,           &
!                                lat, lon, area,                       &
!                                p_half, p_full, z_half, z_full,       &
!                                u, v, t, q, r, um, vm, tm, qm, rm,    &
!                                frac_land, rough_mom,                 &
!                                albedo,    t_surf_rad,                &
!                                u_star,    b_star, q_star,            &
!                                dtau_du,  dtau_dv,  tau_x,  tau_y,    &
!                                udt, vdt, tdt, qdt, rdt,              &
!                                flux_sw,  flux_lw,  coszen,  gust,    &
!                                Surf_diff, gavg_rrv,                  &
!                                mask, kbot
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <INOUT NAME="rd" TYPE="real">
!   multiple 3d diagnostic tracer fields 
!  </INOUT>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <IN NAME="rough_mom" TYPE="real">
!   boundary layer roughness
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="t_surf_rad" TYPE="real">
!   surface radiative temperature
!  </IN>
!  <IN NAME="u_star" TYPE="real">
!   boundary layer wind speed (frictional speed)
!  </IN>
!  <IN NAME="b_star" TYPE="real">
!   ???
!  </IN>
!  <IN NAME="q_star" TYPE="real">
!   boundary layer specific humidity
!  </IN>
!  <IN NAME="dtau_du" TYPE="real">
!   derivative of zonal surface stress w.r.t zonal wind speed
!  </IN>
!  <IN NAME="dtau_dv" TYPE="real">
!   derivative of meridional surface stress w.r.t meridional wind speed
!  </IN>
!  <INOUT NAME="tau_x" TYPE="real">
!   boundary layer meridional component of wind shear
!  </INOUT>
!  <INOUT NAME="tau_y" TYPE="real">
!   boundary layer zonal component of wind shear
!  </INOUT>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="flux_sw" TYPE="real">
!   Shortwave flux from radiation package
!  </OUT>
!  <OUT NAME="flux_lw" TYPE="real">
!   Longwave flux from radiation package
!  </OUT>
!  <OUT NAME="coszen" TYPE="real">
!   cosine of zenith angle
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
!  <IN NAME="gavg_rrv" TYPE="real">
!   array containing global average of tracer volume mixing ratio
!  </IN>
!!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
!
!  <IN NAME="diff_cum_mom" TYPE="real">
!   OPTIONAL: present when do_moist_processes=.false.
!    cu_mo_trans diffusion coefficients, which are passed through to vert_diff_down.
!    Should not be present when do_moist_processes=.true., since these
!    values are passed out from moist_processes.
!  </IN>
!
!  <IN NAME="moist_convect" TYPE="real">
!   OPTIONAL: present when do_moist_processes=.false.
!    Should not be present when do_moist_processes=.true., since these
!    values are passed out from moist_processes.
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_down (is, ie, js, je,                       &
                                Time_prev, Time, Time_next,           &
                                lat, lon, area,                       &
                                p_half, p_full, z_half, z_full,       &
                                u, v, t, q, r, um, vm, tm, qm, rm,    &
                                frac_land, rough_mom,                 &
                                frac_open_sea,                        &
                                albedo, albedo_vis_dir, albedo_nir_dir,&
                                albedo_vis_dif, albedo_nir_dif,       &
                                t_surf_rad,                           &
                                u_star,    b_star, q_star,            &
                                dtau_du, dtau_dv,  tau_x,  tau_y,     &
                                udt, vdt, tdt, qdt, rdt,              &
                                flux_sw,                              &
                                flux_sw_dir,                          &
                                flux_sw_dif,                          &
                                flux_sw_down_vis_dir,                 &
                                flux_sw_down_vis_dif,                 &
                                flux_sw_down_total_dir,               &
                                flux_sw_down_total_dif,               &
                                flux_sw_vis,                          &
                                flux_sw_vis_dir,                      &
                                flux_sw_vis_dif,                      &
                                flux_lw,  coszen,  gust,              &
                                Surf_diff, gavg_rrv,                  &
                                mask, kbot, diff_cum_mom,             &
                                moist_convect, diffm, difft )

!---------------------------------------------------------------------
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!-----------------------------------------------------------------------

integer,                 intent(in)             :: is, ie, js, je
type(time_type),         intent(in)             :: Time_prev, Time,  &
                                                   Time_next
real,dimension(:,:),     intent(in)             :: lat, lon, area
real,dimension(:,:,:),   intent(in)             :: p_half, p_full,   &
                                                   z_half, z_full,   &
                                                   u , v , t , q ,   &
                                                   um, vm, tm, qm
real,dimension(:,:,:,:), intent(inout)          :: r
real,dimension(:,:,:,:), intent(inout)          :: rm
real,dimension(:,:),     intent(in)             :: frac_land,   &
                                                   rough_mom, &
                                                   albedo, t_surf_rad, &
                                                   albedo_vis_dir, albedo_nir_dir, &
                                                   albedo_vis_dif, albedo_nir_dif, &
                                                   u_star, b_star,    &
                                                   q_star, dtau_du,   &
                                                   dtau_dv, frac_open_sea
real,dimension(:,:),     intent(inout)          :: tau_x,  tau_y
real,dimension(:,:,:),   intent(inout)          :: udt,vdt,tdt,qdt
real,dimension(:,:,:,:), intent(inout)          :: rdt
real,dimension(:,:),     intent(out)            :: flux_sw,  &
                                                   flux_sw_dir, &
                                                   flux_sw_dif, flux_lw,  &
                                                   coszen,  gust, &
                                                   flux_sw_down_vis_dir, &
                                                   flux_sw_down_vis_dif, &
                                                   flux_sw_down_total_dir, &
                                                   flux_sw_down_total_dif, &
                                                   flux_sw_vis, &
                                                   flux_sw_vis_dir, & 
                                                   flux_sw_vis_dif 
type(surf_diff_type),    intent(inout)          :: Surf_diff
real,dimension(:),       intent(in)             :: gavg_rrv
real,dimension(:,:,:),   intent(in)   ,optional :: mask
integer, dimension(:,:), intent(in)   ,optional :: kbot
real,  dimension(:,:,:), intent(in)   ,optional :: diff_cum_mom
logical, dimension(:,:), intent(in)   ,optional :: moist_convect
real,  dimension(:,:,:), intent(out)  ,optional :: diffm, difft 
!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      albedo_vis_dir surface visible direct albedo [ dimensionless ]
!      albedo_nir_dir surface nir direct albedo [ dimensionless ]
!      albedo_vis_dif surface visible diffuse albedo [ dimensionless ]
!      albedo_nir_dif surface nir diffuse albedo [ dimensionless ]
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      rd             multiple 3d diagnostic tracer fields 
!                     [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!
!   intent(out) variables:
!
!      flux_sw
!      flux_sw_dir            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_dif            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_down_vis_dir   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_vis_dif   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_total_dir total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_down_total_dif total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_vis            net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dir        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dif        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_lw
!      coszen
!      gust
!
!   intent(in), optional variables:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      real, dimension(size(u,1),size(u,2),size(u,3)) :: diff_t_vert, &
                                                        diff_m_vert, udt_idm, vdt_idm, tdt_idm
      real, dimension(size(rdt,1), size(rdt,2), size(rdt,3), size(rdt,4)) :: rdt_idm
      real, dimension(size(u,1),size(u,2))           :: z_pbl 
      type(aerosol_type)                             :: Aerosol
      type(cld_specification_type)                   :: Cld_spec
      type(radiative_gases_type)                     :: Rad_gases
      type(atmos_input_type)                         :: Atmos_input
      type(surface_type)                             :: Surface
      type(rad_output_type)                          :: Radiation
      type(microphysics_type)                        :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys,&
                                                    Shallow_microphys, &
                                                        Model_microphys
      integer          ::    sec, day, n, nhum
      real             ::    dt, alpha, dt2
      logical          ::    used
      real, dimension(size(u,1),size(u,2))           :: net_surf_sw_down_grey, surf_lw_down_grey, coszen_idmp

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum   
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      Aerosol         aerosol_type variable describing the aerosol
!                      fields to be seen by the radiation package
!      Cld_spec        cld_specification_type variable describing the
!                      cloud field to be seen by the radiation package 
!      Rad_gases       radiative_gases_type variable describing the
!                      radiatively-active gas distribution to be seen 
!                      by the radiation package
!      Atmos_input     atmos_input_type variable describing the atmos-
!                      pheric state to be seen by the radiation package
!      Surface         surface_type variable describing the surface
!                      characteristics to be seen by the radiation 
!                      package
!      Radiation       rad_output_type variable containing the variables
!                      output from the radiation package, for passage
!                      to other modules
!      Rad_time        time at which the radiation calculation is to
!                      apply [ time_type ]
!      Lsc_microphys   microphysics_type variable containing the micro-
!                      physical characteristics of the large-scale
!                      clouds to be seen by the radiation package 
!      Meso_microphys  microphysics_type variable containing the micro-
!                      physical characteristics of the mesoscale
!                      clouds to be seen by the radiation package 
!      Cell_microphys  microphysics_type variable containing the micro-
!                      physical characteristics of the cell-scale
!                      clouds to be seen by the radiation package 
!      Shallow_microphys  
!                      microphysics_type variable containing the micro-
!                      physical characteristics of the cell-scale
!                      clouds to be seen by the radiation package
!      sec, day        second and day components of the time_type 
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!      need_aerosols   need to obtain aerosol data on this time step
!                      to input to the radiation package ?
!      need_clouds     need to obtain cloud data on this time step
!                      to input to the radiation package ?
!      need_gases      need to obtain radiative gas data on this time 
!                      step to input to the radiation package ?
!      need_basic      need to obtain atmospheric state variables on
!                      this time step to input to the radiation package?
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
                         'module has not been initialized', FATAL)
      endif

      if (print_s_messages) write(6,*) 'made it to physics driver down'
!---------------------------------------------------------------------
!    if COSP is activated, save the surface (skin) temperature for
!    its use.
!---------------------------------------------------------------------
      if (do_cosp) then
        tsurf_save(is:ie,js:je) = t_surf_rad
      endif

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
      if (do_check_args) call check_args  &
                   (lat, lon, area, p_half, p_full, z_half, z_full, &
                    u, v, t, q, r, um, vm, tm, qm, rm,              &
                    udt, vdt, tdt, qdt, rdt)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next - Time_prev, sec, day)
      dt = real(sec + day*86400)

!----------------------------------------------------------------------
!    prepare to calculate radiative forcings. obtain the valid time
!    at which the radiation calculation is to apply, the needed atmos-
!    pheric fields, and any needed inputs from other physics modules.
!---------------------------------------------------------------------
      call mpp_clock_begin ( radiation_clock )
      
      
!       PUT RRTM AND GREY RADIATION HERE
      ! subroutine physics_driver_down (is, ie, js, je,                       &
      !   Time_prev, Time, Time_next,           &
      !   lat, lon, area,                       &
      !   p_half, p_full, z_half, z_full,       &
      !   phalfgrey,                            &
      !   u, v, t, q, r, um, vm, tm, qm, rm,    &
      !   frac_land, rough_mom,                 &
      !   frac_open_sea,                        &
      !   albedo, albedo_vis_dir, albedo_nir_dir,&
      !   albedo_vis_dif, albedo_nir_dif,       &
      !   t_surf_rad,                           &
      !   u_star,    b_star, q_star,            &
      !   dtau_du, dtau_dv,  tau_x,  tau_y,     &
      !   udt, vdt, tdt, qdt, rdt,              &
      !   flux_sw,                              &
      !   flux_sw_dir,                          &
      !   flux_sw_dif,                          &
      !   flux_sw_down_vis_dir,                 &
      !   flux_sw_down_vis_dif,                 &
      !   flux_sw_down_total_dir,               &
      !   flux_sw_down_total_dif,               &
      !   flux_sw_vis,                          &
      !   flux_sw_vis_dir,                      &
      !   flux_sw_vis_dif,                      &
      !   flux_lw,  coszen,  gust,              &
      !   Surf_diff, gavg_rrv,                  &
      !   mask, kbot, diff_cum_mom,             &
      !   moist_convect, diffm, difft  )      

      if (print_s_messages) write(6,*) 'maxval tdt pre idm', maxval(tdt), maxval(qdt)

      if (print_s_messages) write(6,*) 'made it to idmp call'
    if (do_grey_radiation) then
      call idealized_radiation_and_optional_surface_flux(is, js, Time, dt, p_half, p_full, z_half, z_full, u, v, t, r, t_surf_rad, udt, vdt, tdt, rdt, .false., net_surf_sw_down_grey=net_surf_sw_down_grey, surf_lw_down_grey = surf_lw_down_grey, coszen_out = coszen_idmp, albedo_in = albedo )  
    else
      call idealized_radiation_and_optional_surface_flux(is, js, Time, dt, p_half, p_full, z_half, z_full, u, v, t, r, t_surf_rad, udt, vdt, tdt, rdt, .false., coszen_out = coszen_idmp, albedo_in = albedo )
    endif
    
!     call get_nhum(nhum)
!     qdt = rdt(:,:,:,nhum)

    if (print_s_messages) write(6,*) 'maxval tdt post idm', maxval(tdt), maxval(qdt)

    if (print_s_messages) write(6,*) 'Done idmp call'
!-------------------------------------------------------------------
!    process the variables returned from radiation_driver_mod. the 
!    radiative heating rate is added to the accumulated physics heating
!    rate (tdt). net surface lw and sw fluxes and the cosine of the 
!    zenith angle are placed in locations where they can be exported
!    for use in other component models. the lw heating rate is stored
!    in a module variable for potential use in other physics modules.
!    the radiative heating rate is also added to a variable which is
!    accumulating the radiative and turbulent heating rates, and which
!    is needed by strat_cloud_mod.
!-------------------------------------------------------------------
    if (print_s_messages) write(6,*) 'applying tendencies'
      ! tdt     = tdt + tdt_idm
    !PROBLEM IS CURRENTLY HERE, where Radiation object is not defined and so getting seg faults. Need to work out relationship between all these variables in order to divide them up. 

    ! real(kind=rb) :: swnflx(nlay+2)         ! Total sky shortwave net flux (W/m2) #NEED THIS
    ! real(kind=rb) :: swnflxc(nlay+2)        ! Clear sky shortwave net flux (W/m2) #DON'T need this
    ! real(kind=rb) :: dirdflux(nlay+2)       ! Direct downward shortwave surface flux #NEED THIS
    ! real(kind=rb) :: difdflux(nlay+2)       ! Diffuse downward shortwave surface flux #DON'T need this
    ! real(kind=rb) :: uvdflx(nlay+2)         ! Total sky downward shortwave flux, UV/vis  #NEED UPWARD AND THIS
    ! real(kind=rb) :: nidflx(nlay+2)         ! Total sky downward shortwave flux, near-IR  #DON'T need this
    ! real(kind=rb) :: dirdnuv(nlay+2)        ! Direct downward shortwave flux, UV/vis #NEED THIS
    ! real(kind=rb) :: difdnuv(nlay+2)        ! Diffuse downward shortwave flux, UV/vis #NEED UPWARD COMPONENT AND THIS
    ! real(kind=rb) :: dirdnir(nlay+2)        ! Direct downward shortwave flux, near-IR #DON'T need this
    ! real(kind=rb) :: difdnir(nlay+2)        ! Diffuse downward shortwave flux, near-IR #NEED UPWARD COMPONENT AND THIS


      ! flux_sw = Radiation%flux_sw_surf(:,:,1) !radiation driver #6026 - net surface sw flux (dfsw - usfw @ surface) ## use swnflx RRTM ##

      ! flux_sw_dir            = Radiation%flux_sw_surf_dir(:,:,1) !radiation driver - dfsw_dir_sfc  ## use dirdflux RRTM    ##

      ! flux_sw_dif            = Radiation%flux_sw_surf_dif(:,:,1) ! dfsw_dif_sfc - ufsw_dif_sfc ## difdnuv+difdnir (total diffuse fluxes)    ##

      ! flux_sw_down_vis_dir   = Radiation%flux_sw_down_vis_dir(:,:,1) ! dfsw_vis_sfc_dir ##  dirdnuv    ##

      ! flux_sw_down_vis_dif   = Radiation%flux_sw_down_vis_dif(:,:,1) !dfsw_vis_sfc_dif ##   difdnuv  ##

      ! flux_sw_down_total_dir = Radiation%flux_sw_down_total_dir(:,:,1) ! dfsw_dir_sfc Seems to be the same as flux_sw_dir ##  dirdflux   ##

      ! flux_sw_down_total_dif = Radiation%flux_sw_down_total_dif(:,:,1) ! dfsw_dif_sfc Just the down component of flux_sw_dif ##  difdnuv+difdnir **NEED TO GET JUST DOWN COMPONENT OF THIS TO MAKE IT WORK**   ##

      ! flux_sw_vis            = Radiation%flux_sw_vis(:,:,1) ! dfsw_vis_sfc - ufsw_vis_sfc ##  uvdflx ** need upward flux component too   ##


      ! flux_sw_vis_dir        = Radiation%flux_sw_vis_dir(:,:,1) ! dfsw_vis_sfc_dir ##  dirdnuv    ##

      
      ! flux_sw_vis_dif        = Radiation%flux_sw_vis_dif(:,:,1) !dfsw_vis_sfc_dif - ufsw_vis_sfc_dif ##  difdnuv  ** NEED UPWARD COMPONENT   ##


      ! flux_lw = Radiation%flux_lw_surf !STEFAN*Atmos_input%temp(:,:,kmax+1)**4 - Lw_output(1)%flxnet(:,:,kmax+1) ##     ##

      ! coszen  = Radiation%coszen_angle !coszen_angle
      ! lw_tendency(is:ie,js:je,:) = Radiation%tdtlw(:,:,:)


      ! radturbten (is:ie,js:je,:) = radturbten(is:ie,js:je,:) + &
      !                              Radiation%tdt_rad(:,:,:,1)
 
      if (print_s_messages) write(6,*) 'DONE applying tendencies'
      call mpp_clock_end ( radiation_clock )

      if (print_s_messages) write(6,*) 'do grey radiation =', do_grey_radiation, maxval(net_surf_sw_down_grey)

      if(do_grey_radiation) then !rif:(09/10/09) 
        ! call grey_radiation(is, js, Time, Time_next, lat, lon, phalfgrey, albedo, t_surf_rad, t, tdt, net_surf_sw_down_grey, flux_lw)
        coszen = coszen_idmp
        flux_sw         = net_surf_sw_down_grey
        flux_sw_dir     = R1*net_surf_sw_down_grey
        flux_sw_dif     = R2*net_surf_sw_down_grey
        flux_sw_vis_dir = R3*net_surf_sw_down_grey
        flux_sw_vis_dif = R4*net_surf_sw_down_grey
        flux_lw         = surf_lw_down_grey
      endif

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired. 
!----------------------------------------------------------------------
      if (print_s_messages) write(6,*) 'about to do damping driver'
      z_pbl(:,:) = pbltop(is:ie,js:je) 
      if (print_s_messages) write(6,*) maxval(pbltop(is:ie, js:je)), 'max pbltop'
      call mpp_clock_begin ( damping_clock )
      call damping_driver (is, js, lat, Time_next, dt,           &
                           p_full, p_half, z_full, z_half,          &
                           um, vm, tm, qm, rm(:,:,:,1:ntp), &
                           udt, vdt, tdt, qdt, rdt,&
                           z_pbl)
      if (print_s_messages) write(6,*) 'done damping driver'
      call mpp_clock_end ( damping_clock )

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------

      ! if (id_tdt_phys_turb > 0) then
      !   used = send_data ( id_tdt_phys_turb, -2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_turb(n) > 0) then
      !     used = send_data ( id_tracer_phys_turb(n), -2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

      if (print_s_messages) write(6,*) 'about to do vert turb driver'
      call mpp_clock_begin ( turb_clock )
      call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             lw_tendency(is:ie,js:je,:), frac_land,  &
                             p_half, p_full, z_half, z_full, u_star, &
                             b_star, q_star, rough_mom, lat,         &
                             convect(is:ie,js:je),                   &
                             u, v, t, q, r(:,:,:,1:ntp), um, vm,     &
                             tm, qm, rm(:,:,:,1:ntp),                &
                             udt, vdt, tdt, qdt, rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl)
     call mpp_clock_end ( turb_clock )
     pbltop(is:ie,js:je) = z_pbl(:,:)
     if (print_s_messages) write(6,*) 'Done vert turb driver'
!      write(6,*) maxval(diff_t_vert), minval(diff_t_vert), maxval(diff_m_vert), minval(diff_m_vert), maxval(z_pbl), minval(z_pbl)
      ! if (id_tdt_phys_turb > 0) then
      !   used = send_data ( id_tdt_phys_turb, +2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_turb(n) > 0) then
      !     used = send_data ( id_tracer_phys_turb(n), +2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------
      ! call mpp_clock_begin ( tracer_clock )
      ! call atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
      !                           area, z_pbl, rough_mom,         &
      !                           frac_open_sea,   &
      !                           frac_land, p_half, p_full,  &
      !                           u, v, t, q, r, &
      !                           rm, rdt, dt, &
      !                           u_star, b_star, q_star, &
      !                           z_half, z_full, t_surf_rad, albedo, &
      !                           Time_next, &
      !                           flux_sw_down_vis_dir, flux_sw_down_vis_dif, &  
      !                           mask, kbot)
      ! call mpp_clock_end ( tracer_clock )

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion 
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current 
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
      if (diffusion_smooth) then
        call get_time (Time_next - Time, sec, day)
        dt2 = real(sec + day*86400)
        alpha = dt2/tau_diff
        diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &
                                 alpha*(diff_m_vert(:,:,:) +  &
                                 diff_cu_mo(is:ie,js:je,:)) )/&
                                 (1. + alpha)
        where (diff_m(is:ie,js:je,:) < diff_min)
          diff_m(is:ie,js:je,:) = 0.0
        end where
        diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                 alpha*diff_t_vert(:,:,:) )/  &
                                 (1. + alpha)
        where (diff_t(is:ie,js:je,:) < diff_min)
          diff_t(is:ie,js:je,:) = 0.0
        end where
      else
        diff_t(is:ie,js:je,:) = diff_t_vert
        diff_m(is:ie,js:je,:) = diff_m_vert + diff_cu_mo(is:ie, js:je,:)
      end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------

      ! if (id_tdt_phys_vdif_dn > 0) then
      !   used = send_data ( id_tdt_phys_vdif_dn, -2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_vdif_dn(n) > 0) then
      !     used = send_data ( id_tracer_phys_vdif_dn(n), -2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

      call mpp_clock_begin ( diff_down_clock )
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:)
    !   if ( allocated(diff_t_clubb) ) then  !RASF Bug fix. Can only pass bounds of diff_t_clubb when allocated
    !      call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
    !                               p_full, z_full,   &
    !                               diff_m(is:ie,js:je,:),         &
    !                               diff_t(is:ie,js:je,:),         &
    !                               um ,vm ,tm ,qm ,rm(:,:,:,1:ntp), &
    !                               dtau_du, dtau_dv, tau_x, tau_y,  &
    !                               udt, vdt, tdt, qdt, rdt,       &
    !                               Surf_diff,                     &
    !                               diff_t_clubb=diff_t_clubb(is:ie,js:je,:),   &   ! cjg
    !                               mask=mask, kbot=kbot           )
    !  else
      ! write(6,*) 'max in tau_x, tau_y', maxval(tau_x), minval(tau_x), maxval(tau_y), minval(tau_y)
      ! write(6,*) 'max in dt_tg', maxval(tdt), minval(tdt)

        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                  p_full, z_full,   &
                                  diff_m(is:ie,js:je,:),         &
                                  diff_t(is:ie,js:je,:),         &
                                  um ,vm ,tm ,qm ,rm(:,:,:,1:ntp), &
                                  dtau_du, dtau_dv, tau_x, tau_y,  &
                                  udt, vdt, tdt, qdt, rdt,       &
                                  Surf_diff           )
    !   endif
      ! write(6,*) 'max out dt_tg', maxval(tdt), minval(tdt)
      ! write(6,*) 'max out tau_x, tau_y', maxval(tau_x), minval(tau_x), maxval(tau_y), minval(tau_y)

      ! if (id_tdt_phys_vdif_dn > 0) then
      !   used = send_data ( id_tdt_phys_vdif_dn, +2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_vdif_dn(n) > 0) then
      !     used = send_data ( id_tracer_phys_vdif_dn(n), +2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t(is:ie,js:je,:)
      endif
      if (present(diffm)) then
        diffm = diff_m(is:ie,js:je,:)
      endif

     call mpp_clock_end ( diff_down_clock )

     if (print_s_messages) write(6,*) 'maxval tdt end pdd', maxval(tdt), maxval(qdt)

     if (print_s_messages) write(6,*) 'Finished physics driver down'
 end subroutine physics_driver_down



!#######################################################################
! <SUBROUTINE NAME="physics_driver_up">
!  <OVERVIEW>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_up (is, ie, js, je,                    &
!                               Time_prev, Time, Time_next,        &
!                               lat, lon, area,                    &
!                               p_half, p_full, z_half, z_full,    & 
!                               omega,                             &
!                               u, v, t, q, r, um, vm, tm, qm, rm, &
!                               frac_land,                         &
!                               udt, vdt, tdt, qdt, rdt,           &
!                               Surf_diff,                         &
!                               lprec,   fprec, gust,              &
!                               mask, kbot    )
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="omega" TYPE="real">
!   Veritical pressure tendency
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="lprec" TYPE="real">
!  </OUT>
!  <OUT NAME="fprec" TYPE="real">
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
 subroutine physics_driver_up (is, ie, js, je,                    &
                               Time_prev, Time, Time_next,        &
                               lat, lon, area,                    &
                               p_half, p_full, z_half, z_full,    &
                               omega,                             &
                               u, v, t, q, r, um, vm, tm, qm, rm, &
                               frac_land,                         &
                               u_star, b_star, q_star,            &
                               udt, vdt, tdt, qdt, rdt,           &
                               Surf_diff,                         &
                               lprec,   fprec, gust,              &
                               mask, kbot,                        &
                               hydrostatic, phys_hydrostatic)

!----------------------------------------------------------------------
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!---------------------------------------------------------------------

integer,                intent(in)             :: is, ie, js, je
type(time_type),        intent(in)             :: Time_prev, Time,   &
                                                  Time_next
real,dimension(:,:),    intent(in)             :: lat, lon, area
real,dimension(:,:,:),  intent(in)             :: p_half, p_full,   &
                                                  omega,  &
                                                  z_half, z_full,     &
                                                  u , v , t , q ,    &
                                                  um, vm, tm, qm
real,dimension(:,:,:,:),intent(in)          :: r,rm       ! cjg: inout
real,dimension(:,:),    intent(in)             :: frac_land
real,dimension(:,:),    intent(in)             :: u_star, b_star, q_star
real,dimension(:,:,:),  intent(inout)          :: udt,vdt,tdt,qdt
real,dimension(:,:,:,:),intent(inout)          :: rdt
type(surf_diff_type),   intent(inout)          :: Surf_diff
real,dimension(:,:),    intent(out)            :: lprec, fprec
real,dimension(:,:),    intent(inout)          :: gust
real,dimension(:,:,:),  intent(in),   optional :: mask
integer,dimension(:,:), intent(in),   optional :: kbot
logical,                intent(in),   optional :: hydrostatic, phys_hydrostatic

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      omega
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!      gust
!
!   intent(out) variables:
!
!      lprec     
!      fprec       
!
!   intent(in), optional variables:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!   local variables:

      real, dimension(size(u,1), size(u,2), size(u,3)) :: diff_cu_mo_loc
      real, dimension(size(u,1), size(u,2))            :: gust_cv
      real, dimension(size(u,1), size(u,2))            :: land_mask
      integer :: sec, day
      real    :: dt
      real, dimension(size(t,1), size(t,2)) :: u_sfc, v_sfc
      real, dimension(size(t,1), size(t,2), size(t,3)+1) :: pflux
      real, dimension(size(t,1), size(t,2), size(t,3))   ::  &
                             tca, cca, rhoi, lsliq, lsice, ccliq,  &
                             ccice, reff_lsclliq, reff_lsclice, &
                             reff_ccclliq, reff_ccclice, &
                             reff_lsprliq, reff_lsprice, &
                             reff_ccprliq, reff_ccprice, &
                             fl_lsrain_loc, fl_lssnow_loc,  &
                             fl_lsgrpl_loc, &
                             fl_donmca_rain_loc, fl_donmca_snow_loc, &
                             fl_ccrain_loc, fl_ccsnow_loc, mr_ozone_loc
      real, dimension(size(t,1), size(t,2), size(t,3), ncol) ::  &
                             stoch_mr_liq, stoch_mr_ice, &
                             stoch_size_liq, stoch_size_frz
      integer :: i, j , k, n
      integer :: nls, ncc
      real    :: alphb
      integer :: flag_ls, flag_cc
      integer :: kmax
      logical :: used

! ---> h1g, 2012-08-28,  
! save the temperature and moisture tendencies from sensible and latent heat fluxes
  real, dimension(size(t,1), size(t,2))  ::   tdt_shf,  qdt_lhf
! <--- h1g, 2012-08-28

!---------------------------------------------------------------------
!   local variables:
!
!        diff_cu_mo_loc   diffusion coefficient contribution due to 
!                         cumulus momentum transport
!        gust_cv
!        sec, day         second and day components of the time_type 
!                         variable
!        dt               physics time step [ seconds ]
!
!---------------------------------------------------------------------
      ! type(aerosol_type)                               :: Aerosol

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
             'module has not been initialized', FATAL)
      endif

!----------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kmax = size(u,3)

!----------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------

      ! if (id_tdt_phys_vdif_up > 0) then
      !   used = send_data ( id_tdt_phys_vdif_up, -2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_vdif_up(n) > 0) then
      !     used = send_data ( id_tracer_phys_vdif_up(n), -2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

      call mpp_clock_begin ( diff_up_clock )
      call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                Surf_diff, tdt, qdt, rdt)

      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt(:,:,:)
      call mpp_clock_end ( diff_up_clock )

      ! if (id_tdt_phys_vdif_up > 0) then
      !   used = send_data ( id_tdt_phys_vdif_up, +2.0*tdt(:,:,:), &
      !                      Time_next, is, js, 1, rmask=mask )
      ! endif

      ! do n=1,nt
      !   if (id_tracer_phys_vdif_up(n) > 0) then
      !     used = send_data ( id_tracer_phys_vdif_up(n), +2.0*rdt(:,:,:,n), &
      !                        Time_next, is, js, 1, rmask=mask )
      !   endif
      ! end do

!-----------------------------------------------------------------------
!    if the fms integration path is being followed, call moist processes
!    to compute moist physics, including convection and processes 
!    involving condenstion.
!-----------------------------------------------------------------------
      if (do_moist_processes) then

        ! if (id_tdt_phys_moist > 0) then
        !   used = send_data ( id_tdt_phys_moist, -2.0*tdt(:,:,:), &
        !                      Time_next, is, js, 1, rmask=mask )
        ! endif

        ! do n=1,nt
        !   if (id_tracer_phys_moist(n) > 0) then
        !     used = send_data ( id_tracer_phys_moist(n), -2.0*rdt(:,:,:,n), &
        !                        Time_next, is, js, 1, rmask=mask )
        !   endif
        ! end do

        call mpp_clock_begin ( moist_processes_clock )
        if (print_s_messages) write(6,*) 'about to do moist processes'
        call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                           p_half, p_full, z_half, z_full, omega,    &
                           diff_t(is:ie,js:je,:),                    &
                           radturbten(is:ie,js:je,:),                &
                           cush           (is:ie,js:je),          &!
                        cbmf           (is:ie,js:je),             &!
                            pbltop(is:ie,js:je),         &!miz
                            u_star, b_star, q_star,          &!miz
                           t, q, r, u, v, tm, qm, rm, um, vm,        &
                           tdt, qdt, rdt, udt, vdt, diff_cu_mo_loc , &
                           convect(is:ie,js:je), lprec, fprec,       &
               fl_lsrain(is:ie,js:je,:), fl_lssnow(is:ie,js:je,:),  &
               fl_ccrain(is:ie,js:je,:), fl_ccsnow(is:ie,js:je,:),    &
           fl_donmca_rain(is:ie,js:je,:), fl_donmca_snow(is:ie,js:je,:), &
                           gust_cv, area, lon, lat)      
                           
        if (print_s_messages) write(6,*) 'Done moist processes'                           
        call mpp_clock_end ( moist_processes_clock )
        diff_cu_mo(is:ie, js:je,:) = diff_cu_mo_loc(:,:,:)
        radturbten(is:ie,js:je,:) = 0.0

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained 
!    from non-convective parameterizations.
!---------------------------------------------------------------------
        gust = sqrt( gust*gust + gust_cv*gust_cv)

!         if (id_tdt_phys_moist > 0) then
!           used = send_data ( id_tdt_phys_moist, +2.0*tdt(:,:,:), &
!                              Time_next, is, js, 1, rmask=mask )
!         endif

!         do n=1,nt
!           if (id_tracer_phys_moist(n) > 0) then
!             used = send_data ( id_tracer_phys_moist(n), +2.0*rdt(:,:,:,n), &
!                                Time_next, is, js, 1, rmask=mask )
!           endif
!         end do

!         if (id_tdt_phys > 0) then
!            used = send_data ( id_tdt_phys, tdt(:,:,:), &
!                               Time_next, is, js, 1, rmask=mask )
!         endif
! !--> cjg
! !       if (id_qdt_phys > 0) then
! !          used = send_data ( id_qdt_phys, qdt(:,:,:), &
! !                             Time_next, is, js, 1, rmask=mask )
! !       endif
!         do n=1,nt
!           if (id_tracer_phys(n) > 0) then
!             used = send_data ( id_tracer_phys(n), rdt(:,:,:,n), &
!                                Time_next, is, js, 1, rmask=mask )
!           endif
!         end do
!<--cjg
      endif


 end subroutine physics_driver_up


!#######################################################################
! <SUBROUTINE NAME="physics_driver_end">
!  <OVERVIEW>
!   physics_driver_end is the destructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_end is the destructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_end (Time)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_end (Time)

!---------------------------------------------------------------------
!    physics_driver_end is the destructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      Time      current time [ time_type(days, seconds) ]
!
!--------------------------------------------------------------------
integer :: moist_processes_term_clock, damping_term_clock, turb_term_clock, &
           diff_term_clock, cloud_spec_term_clock, aerosol_term_clock, &
           grey_radiation_term_clock, radiative_gases_term_clock, &
           radiation_term_clock, tracer_term_clock, cosp_term_clock

! ---> h1g
integer :: clubb_term_clock

write(6,*) 'physics driver end #0'

!--------------------------------------------------------------------
!         clubb_term_clock =      &
!         mpp_clock_id( '   Phys_driver_term: clubb: Termination', &
!                 grain=CLOCK_MODULE_DRIVER )
! ! <--- h1g

!       moist_processes_term_clock =      &
!         mpp_clock_id( '   Phys_driver_term: MP: Termination', &
!                 grain=CLOCK_MODULE_DRIVER )
!       damping_term_clock         =     &
!         mpp_clock_id( '   Phys_driver_term: Damping: Termination',    &
!                   grain=CLOCK_MODULE_DRIVER )
!       turb_term_clock            =      &
!         mpp_clock_id( '   Phys_driver_term: Vert. Turb.: Termination', &
!                   grain=CLOCK_MODULE_DRIVER )
!       diff_term_clock       =     &
!         mpp_clock_id( '   Phys_driver_term: Vert. Diff.: Termination',   &
!                  grain=CLOCK_MODULE_DRIVER )
!       cloud_spec_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: Cloud spec: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       cosp_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: COSP: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       aerosol_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: Aerosol: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       grey_radiation_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: Grey Radiation: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       radiative_gases_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: Radiative gases: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       radiation_term_clock       =       &
!         mpp_clock_id( '   Phys_driver_term: Radiation: Termination', &
!                        grain=CLOCK_MODULE_DRIVER )
!       tracer_term_clock          =      &
!         mpp_clock_id( '   Phys_driver_term: Tracer: Termination',    &
!                  grain=CLOCK_MODULE_DRIVER )
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
              'module has not been initialized', FATAL)
      endif

      write(6,*) 'physics driver end #1'


      ! call physics_driver_netcdf

!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
      ! call mpp_clock_begin ( turb_term_clock )
      call vert_turb_driver_end
      ! call mpp_clock_end ( turb_term_clock )
      write(6,*) 'physics driver end #1.1'


      ! call mpp_clock_begin ( diff_term_clock )
      ! call vert_diff_driver_end
      ! call mpp_clock_end ( diff_term_clock )
      if (do_radiation) then
        ! call mpp_clock_begin ( radiation_term_clock )
        ! call radiation_driver_end
        ! call mpp_clock_end ( radiation_term_clock )
        ! call mpp_clock_begin ( radiative_gases_term_clock )
        ! call radiative_gases_end
        write(6,*) 'physics driver end #1.2'

        ! call mpp_clock_end ( radiative_gases_term_clock )
        ! call mpp_clock_begin ( cloud_spec_term_clock )
        ! call cloud_spec_end
        ! call mpp_clock_end ( cloud_spec_term_clock )
        ! call mpp_clock_begin ( aerosol_term_clock )
        ! call aerosol_end
        ! call mpp_clock_end ( aerosol_term_clock )
      endif
      write(6,*) 'physics driver end #1.3'

      ! call mpp_clock_begin ( grey_radiation_term_clock )

      write(6,*) 'physics driver end #2'

      ! if(do_grey_radiation) call grey_radiation_end 

      ! call mpp_clock_end ( grey_radiation_term_clock )
      ! call mpp_clock_begin ( moist_processes_term_clock )
      call moist_processes_end( clubb_term_clock )
      ! call mpp_clock_end ( moist_processes_term_clock )
      ! call mpp_clock_begin ( tracer_term_clock )
      ! call atmos_tracer_driver_end
      ! call mpp_clock_end ( tracer_term_clock )
      ! call mpp_clock_begin ( damping_term_clock )
      call damping_driver_end
      ! call mpp_clock_end ( damping_term_clock )
      ! call mpp_clock_begin ( cosp_term_clock )
      ! call mpp_clock_end ( cosp_term_clock )

      write(6,*) 'physics driver end #3'      
!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
      deallocate (diff_cu_mo, diff_t, diff_m, pbltop, cush, cbmf,  &
                  convect, radturbten, lw_tendency, r_convect)
      deallocate (fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow, &
                  fl_donmca_snow, fl_donmca_rain, mr_ozone, daytime, &
                  temp_last, q_last)
      ! deallocate (lsc_liquid, lsc_ice, &
                  ! lsc_droplet_number, lsc_ice_number)
      ! deallocate (lsc_snow, lsc_rain, lsc_snow_size, lsc_rain_size)
      if (do_clubb > 0) then
         deallocate ( dcond_ls_liquid )
         deallocate ( dcond_ls_ice )
         deallocate ( Ndrop_act_CLUBB )
         deallocate ( Icedrop_act_CLUBB )
         deallocate ( ndust )
         deallocate ( rbar_dust )
         deallocate ( diff_t_clubb )
      end if

      write(6,*) 'physics driver end #4'

      if (doing_donner) then
        deallocate (cell_cld_frac, cell_liq_amt, cell_liq_size, &
                    cell_ice_amt, cell_ice_size, cell_droplet_number, &
                    meso_cld_frac, meso_liq_amt, meso_liq_size, &
                    meso_ice_amt, meso_ice_size, meso_droplet_number, &
                    nsum_out)
      endif
      if (doing_uw_conv) then
        deallocate (shallow_cloud_area, shallow_liquid, shallow_ice, &
                    shallow_droplet_number)
        deallocate ( shallow_ice_number)

      endif
 
      ! deallocate (id_tracer_phys_vdif_dn)
      ! deallocate (id_tracer_phys_vdif_up)
      ! deallocate (id_tracer_phys_turb)
      ! deallocate (id_tracer_phys_moist)

      write(6,*) 'physics driver end #5'


      ! if (do_cosp .or. do_modis_yim) then
      !   deallocate (stoch_cloud_type, tau_stoch, lwem_stoch, &
      !               stoch_conc_drop, stoch_conc_ice, stoch_size_drop, &
      !               stoch_size_ice, tsurf_save)
      ! endif
!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!-----------------------------------------------------------------------

 end subroutine physics_driver_end

!#######################################################################
! <SUBROUTINE NAME="physics_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp


  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
  endif
  call physics_driver_netcdf(timestamp)
  call vert_turb_driver_restart(timestamp)
  if (do_radiation) then
    ! call radiation_driver_restart(timestamp)
    ! call radiative_gases_restart(timestamp)
  endif

!    call moist_processes_restart(timestamp)
  ! call damping_driver_restart(timestamp)

end subroutine physics_driver_restart
! </SUBROUTINE> NAME="physics_driver_restart"

! <SUBROUTINE NAME="physics_driver_netcdf">
!
! <DESCRIPTION>
! Write out restart file for physics driver.
! This routine is needed so that physics_driver_restart and physics_driver_end
! can call a routine which will not result in multiple copies of restart files 
! being written by the destructor routines.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_netcdf(timestamp)
  character(len=*), intent(in), optional :: timestamp

    r_convect = 0.
    where(convect)
       r_convect = 1.0
    end where
    call save_restart(Phy_restart, timestamp)
    if(in_different_file) call save_restart(Til_restart, timestamp)

end subroutine physics_driver_netcdf
! </SUBROUTINE> NAME="physics_driver_netcdf"

!#######################################################################
! <FUNCTION NAME="do_moist_in_phys_up">
!  <OVERVIEW>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </OVERVIEW>
!  <DESCRIPTION>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </DESCRIPTION>
!  <TEMPLATE>
!   logical = do_moist_in_phys_up()
!  </TEMPLATE>
! </FUNCTION>
!
function do_moist_in_phys_up()

!--------------------------------------------------------------------
!    do_moist_in_phys_up returns the value of do_moist_processes
!----------------------------------------------------------------------

logical :: do_moist_in_phys_up

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('do_moist_in_phys_up',  &
              'module has not been initialized', FATAL)
      endif
 
!-------------------------------------------------------------------
!    define output variable.
!-------------------------------------------------------------------
      do_moist_in_phys_up = do_moist_processes

 
end function do_moist_in_phys_up

!#####################################################################
! <FUNCTION NAME="get_diff_t">
!  <OVERVIEW>
!    returns the values of array diff_t
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array diff_t
!  </DESCRIPTION>
!  <TEMPLATE>
!   diff_t(:,:,:) = get_diff_t()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_diff_t() result(diff_t_out)
real, dimension(size(diff_t,1),size(diff_t,2),size(diff_t,3)) :: diff_t_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_diff_t','module has not been initialized', FATAL)
  endif

  diff_t_out = diff_t

end function get_diff_t

!#####################################################################
! <FUNCTION NAME="get_radturbten">
!  <OVERVIEW>
!    returns the values of array radturbten
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array radturbten
!  </DESCRIPTION>
!  <TEMPLATE>
!   radturbten(:,:,:) = get_radturbten()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_radturbten() result(radturbten_out)
real, dimension(size(radturbten,1),size(radturbten,2),size(radturbten,3)) :: radturbten_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_radturbten','module has not been initialized', FATAL)
  endif

  radturbten_out = radturbten

end function get_radturbten
!#####################################################################
! <SUBROUTINE NAME="zero_radturbten">
!  <OVERVIEW>
!    sets all values of array radturbten to zero
!  </OVERVIEW>
!  <DESCRIPTION>
!    sets all values of array radturbten to zero
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zero_radturbten()
!  </TEMPLATE>
! </SUBROUTINE>
!
!#####################################################################
subroutine zero_radturbten()

  if ( .not. module_is_initialized) then
    call error_mesg ('zero_radturbten','module has not been initialized', FATAL)
  endif

  radturbten = 0.0

end subroutine zero_radturbten
!#####################################################################



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
               
     
!#####################################################################
! <SUBROUTINE NAME="physics_driver_register_restart">
!  <OVERVIEW>
!    physics_driver_register_restart will register restart field when do_netcdf file 
!    is true. 
!  </OVERVIEW>
subroutine physics_driver_register_restart
  character(len=64) :: fname, fname2
  integer           :: id_restart

  
  if(doing_strat()) then 
     now_doing_strat = 1
  else
     now_doing_strat = 0
  endif

  if(doing_edt) then 
     now_doing_edt = 1
  else
     now_doing_edt = 0
  endif

  if(doing_entrain) then 
     now_doing_entrain = 1
  else
     now_doing_entrain = 0
  endif

  fname = 'physics_driver.res.nc'
  call get_mosaic_tile_file(fname, fname2, .false. ) 
  allocate(Phy_restart)
  if(trim(fname2) == trim(fname)) then
     Til_restart => Phy_restart
     in_different_file = .false.
  else
     in_different_file = .true.
     allocate(Til_restart)
  endif

  id_restart = register_restart_field(Phy_restart, fname, 'vers',          vers,              no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_strat',   now_doing_strat,   no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_edt',     now_doing_edt,     no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_entrain', now_doing_entrain, no_domain=.true.)

  id_restart = register_restart_field(Til_restart, fname, 'diff_cu_mo', diff_cu_mo)
  id_restart = register_restart_field(Til_restart, fname, 'pbltop',     pbltop)
  id_restart = register_restart_field(Til_restart, fname, 'cush',       cush, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cbmf',       cbmf, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'diff_t',     diff_t)
  id_restart = register_restart_field(Til_restart, fname, 'diff_m',     diff_m)
!-->cjg
  if (do_clubb > 0) then
     id_restart = register_restart_field(Til_restart, fname, 'diff_t_clubb', diff_t_clubb, mandatory = .false.)
  end if
!<--cjg
  id_restart = register_restart_field(Til_restart, fname, 'convect',    r_convect) 
  if (doing_strat()) then
     id_restart = register_restart_field(Til_restart, fname, 'radturbten',       radturbten)
  endif
  if (doing_edt .or. doing_entrain) then
     id_restart = register_restart_field(Til_restart, fname, 'lw_tendency',      lw_tendency)
  endif
  if (doing_donner) then
     id_restart = register_restart_field(Til_restart, fname, 'cell_cloud_frac',  cell_cld_frac, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_amt',  cell_liq_amt,  mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_size', cell_liq_size, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'cell_ice_amt',     cell_ice_amt,  mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'cell_ice_size',    cell_ice_size, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'meso_cloud_frac',  meso_cld_frac, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_amt',  meso_liq_amt,  mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_size', meso_liq_size, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'meso_ice_amt',     meso_ice_amt,  mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'meso_ice_size',    meso_ice_size, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'nsum',             nsum_out,      mandatory = .false.)
  endif
  if (doing_uw_conv) then
     id_restart = register_restart_field(Til_restart, fname, 'shallow_cloud_area',     shallow_cloud_area,     mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'shallow_liquid',         shallow_liquid,         mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'shallow_ice',            shallow_ice,            mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'shallow_droplet_number', shallow_droplet_number, mandatory = .false.)
     id_restart = register_restart_field(Til_restart, fname, 'shallow_ice_number',     shallow_ice_number,     mandatory = .false.)
  endif
  id_restart = register_restart_field(Til_restart, fname, 'lsc_cloud_area',     lsc_cloud_area,     mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_liquid',         lsc_liquid ,        mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_ice',            lsc_ice,            mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_droplet_number', lsc_droplet_number, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_ice_number',     lsc_ice_number,     mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_snow',           lsc_snow,           mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_rain',           lsc_rain,           mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_snow_size',      lsc_snow_size,      mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'lsc_rain_size',      lsc_rain_size,      mandatory = .false.)

end subroutine physics_driver_register_restart
! </SUBROUTINE>    
!#####################################################################
! <SUBROUTINE NAME="check_args">
!  <OVERVIEW>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
!                        u, v, t, q, r, um, vm, tm, qm, rm,             &
!                        udt, vdt, tdt, qdt, rdt, mask, kbot)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="udt" TYPE="real">
!   zonal wind tendency
!  </IN>
!  <IN NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </IN>
!  <IN NAME="tdt" TYPE="real">
!   temperature tendency
!  </IN>
!  <IN NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </IN>
!  <IN NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, mask, kbot)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

real,    dimension(:,:),    intent(in)          :: lat, lon, area
real,    dimension(:,:,:),  intent(in)          :: p_half, p_full,   &
                                                   z_half, z_full,   &
                                                   u, v, t, q, um, vm, &
                                                   tm, qm
real,    dimension(:,:,:,:),intent(in)          :: r, rm
real,    dimension(:,:,:),  intent(in)          :: udt, vdt, tdt, qdt
real,    dimension(:,:,:,:),intent(in)          :: rdt
real,    dimension(:,:,:),  intent(in),optional :: mask
integer, dimension(:,:),    intent(in),optional :: kbot

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer ::  id, jd, kd  ! model dimensions on the processor  
      integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
      id = size(u,1) 
      jd = size(u,2) 
      kd = size(u,3) 

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

      if (nt > 0) then
        ierr = ierr + check_dim (r,  'r',   id,jd,kd,nt)
        ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,nt)
      endif
      if (ntp > 0) then
        ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,ntp)
      endif

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
      if (ierr > 0) then
        call error_mesg ('physics_driver_mod', 'bad dimensions', FATAL)
      endif

!-----------------------------------------------------------------------


      end subroutine check_args


!#######################################################################
! <FUNCTION NAME="check_dim_2d">
!  <OVERVIEW>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_2d (data,name,id,jd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd" TYPE="integer">
!   expected i and j dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_2d (data,name,id,jd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:) :: data
character(len=*), intent(in)        :: name
integer, intent(in)                 :: id, jd
integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
             'dimension 1 of argument ' //  &
              name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

!----------------------------------------------------------------------

      end function check_dim_2d

!#######################################################################
! <FUNCTION NAME="check_dim_3d">
!  <OVERVIEW>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_3d (data,name,id,jd, kd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd" TYPE="integer">
!   expected i, j and k dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_3d (data,name,id,jd,kd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:,:) :: data
character(len=*), intent(in)          :: name
integer, intent(in)                   :: id, jd, kd
integer  ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
              'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_3d


!#######################################################################
! <FUNCTION NAME="check_dim_4d">
!  <OVERVIEW>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_4d (data,name,id,jd, kd, nt) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd, nt" TYPE="integer">
!   expected i, j, k and 4th dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
real,    intent(in), dimension(:,:,:,:) :: data
character(len=*), intent(in)            :: name
integer, intent(in)                     :: id, jd, kd, nt
integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_4d



!#######################################################################
 

 
                end module physics_driver_mod
