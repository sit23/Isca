
                    module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             relaxed arakawa-schubert
!             donner deep convection
!             large-scale condensation
!             stratiform prognostic cloud scheme 
!             rel humidity cloud scheme 
!             diagnostic cloud scheme 
!             lin cloud microphysics
!             betts-miller convective adjustment
!
!-----------------------------------------------------------------------

! fms modules
use sat_vapor_pres_mod,    only: compute_qs, lookup_es
use time_manager_mod,      only: time_type, get_time
use diag_manager_mod,      only: register_diag_field, send_data
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number,           &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 MPP_CLOCK_SYNC, read_data, write_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_number_tracers, &
                                 get_tracer_names, &
                                 query_method, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 RDGAS, RVGAS, TFREEZE, WTMAIR, &
                                 SECONDS_PER_DAY, KAPPA
! atmos_param modules
! use betts_miller_mod,      only: betts_miller, betts_miller_init
! use bm_massflux_mod,       only: bm_massflux, bm_massflux_init
! use bm_omp_mod,            only: bm_omp, bm_omp_init
! use donner_deep_mod,       only: donner_deep_init,               &
!                                  donner_deep_time_vary,  &
!                                  donner_deep_endts,         &
!                                  donner_deep, donner_deep_end,   &
!                                  donner_deep_restart
! use moist_conv_mod,        only: moist_conv, moist_conv_init
! use lscale_cond_mod,       only: lscale_cond_init
! use uw_conv_mod,           only: uw_conv_end, uw_conv_init
! use lin_cld_microphys_mod, only: lin_cld_microphys_init, &
!                                  lin_cld_microphys_end
! use ras_mod,               only: ras_end, ras_init
! use dry_adj_mod,           only: dry_adj, dry_adj_init
! use strat_cloud_mod,       only: strat_cloud_init, strat_cloud_end, &
!                                  strat_cloud_restart, strat_cloud_time_vary
! use detr_ice_num_mod,      only: detr_ice_num, detr_ice_num_init,   &
!                                  detr_ice_num_end

! ---> h1g
use mpp_mod,               only: mpp_chksum
! use MG_microp_3D_mod,      only: MG_microp_3D_init, MG_microp_3D, &
!                                             MG_microp_3D_end
! 
! use clubb_driver_mod,      only: clubb_init, clubb, clubb_end
! <--- h1g

! use rh_clouds_mod,         only: rh_clouds_init, rh_clouds_end, &
!                                  rh_clouds_sum
! use diag_cloud_mod,        only: diag_cloud_init, diag_cloud_end, &
!                                  diag_cloud_restart
use diag_integral_mod,     only: diag_integral_field_init, &
                                 sum_diag_integral_field
! use cu_mo_trans_mod,       only: cu_mo_trans_init, cu_mo_trans, cu_mo_trans_end
! use moz_hook_mod,          only: moz_hook
use rad_utilities_mod,     only: aerosol_type
use moist_proc_utils_mod,  only: capecalcnew, tempavg, column_diag, rh_calc, pmass

use idealized_moist_phys_mod, only: idealized_convection_and_lscale_cond

! use moistproc_kernels_mod, only: moistproc_init, moistproc_end, moistproc_mca, &
!                                  moistproc_ras, moistproc_lscale_cond,         &
!                                  moistproc_strat_cloud, moistproc_cmt,         &
!                                  moistproc_uw_conv, moistproc_scale_uw,        &
!                                  moistproc_scale_donner,                       &
!                                  rain_uw, snow_uw, ttnd_uw, qtnd_uw, utnd_uw,  &
!                                  vtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw,        &
!                                  qntnd_uw, qtruw, qlin, qiin, qain, delta_ql,  &
!                                  delta_qi, delta_qa, qnitnd_uw
! atmos_shared modules
! use atmos_tracer_utilities_mod, only : wet_deposition

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end, &
            moist_alloc_init, moist_alloc_end,  set_cosp_precip_sources, &
            moist_processes_endts, &
            doing_strat, moist_processes_restart
  

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
   character(len=128) :: &
   version = '$Id: moist_processes.F90,v 20.0 2013/12/13 23:18:25 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal $'

   character(len=5), private :: mod_name = 'moist'
   logical            :: moist_allocated = .false.
   logical            :: module_is_initialized = .false.

!-------------------- namelist data (private) --------------------------

!---------------- namelist variable definitions ------------------------
!
!   do_limit_donner = limit Donner deeo tendencies to prevent the
!                formation of grid points with negative water vapor,
!                liquid or ice.
!
!   do_limit_uw = limit UW shallow tendencies to prevent the formation
!                of grid points with negative total water specific 
!                humidities. This situation can occur because both
!                shallow and deep convection operate on the same
!                soundings without knowledge of what the other is doing
!
!   do_unified_convective_closure = use cloud base mass flux calculated
!                in uw_conv module as value for donner deep parameter-
!                ization; adjust cbmf available for uw shallow appropr-
!                iately. only available when uw shallow and donner deep
!                are the active convective schemes
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=true ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true ]
!   do_ras   = switch to turn on/off relaxed arakawa shubert
!                [logical, default: do_ras=false ]
!   do_donner_deep = switch to turn on/off donner deep convection scheme
!                [logical, default: do_donner_deep=false ]
!   do_strat = switch to turn on/off stratiform cloud scheme
!                [logical, default: do_strat=false ]
!   do_rh_clouds = switch to turn on/off simple relative humidity cloud scheme
!                [logical, default: do_rh_clouds=false ]
!   do_diag_clouds = switch to turn on/off (Gordon's) diagnostic cloud scheme
!                [logical, default: do_diag_clouds=false ]
!   do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=false ]
!   do_lin_cld_microphys = switch to turn on/off the Lin Cloud Micro-Physics scheme
!                [logical, default: do_lin_cld_microphys=false ]
!   do_liq_num = switch to turn on/off the prognostic droplet number scheme.
!                [logical, default: do_liq_num=false ]
!   use_tau  = switch to determine whether current time level (tau)
!                will be used or else future time level (tau+1).
!                if use_tau = true then the input values for t,q, and r
!                are used; if use_tau = false then input values
!                tm+tdt*dt, etc. are used.
!                [logical, default: use_tau=false ]
!
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!
!   do_gust_cv = switch to use convective gustiness (default = false)
!   gustmax    = maximum convective gustiness (m/s)
!   gustconst  = precip rate which defines precip rate which begins to
 !               matter for convective gustiness (kg/m2/sec)
!   cmt_mass_flux_source = parameterization(s) being used to supply the 
!                mass flux profiles seen by the cumulus momentum transport
!                module; currently either 'ras', 'donner', 'uw', 
!                'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 
!                'donner_and_ras_and_uw' or 'all'
!
!   do_bm    = switch to turn on/off betts-miller scheme
!                [logical, default: do_bm=false ]
!   do_bmmass  = switch to turn on/off betts-miller massflux scheme
!                [logical, default: do_bmmass=false ]
!   do_bmomp  = switch to turn on/off olivier's version of the betts-miller 
!                scheme (with separated boundary layer)
!                [logical, default: do_bmomp=false ]
!   do_simple = switch to turn on alternative definition of specific humidity.
!                When true, specific humidity = (rdgas/rvgas)*esat/pressure
!
!   notes: 1) do_lsc and do_strat cannot both be true
!          2) pdepth and tfreeze are used to determine liquid vs. solid
!             precipitation for mca, lsc, and ras schemes, the 
!             stratiform scheme determines it's own precipitation type.
!          3) if do_strat=true then stratiform cloud tracers: liq_wat,
!             ice_wat, cld_amt must be present 
!          4) do_donner_deep and do_rh_clouds cannot both be true
!             (pending revision of code flow)
!
!-----------------------------------------------------------------------
! main convection/large-scale schemes
   logical :: do_bm=.false.
   logical :: do_bmmass =.false.
   logical :: do_bmomp  =.false.
   logical :: do_cmt=.false.
   logical :: do_diag_clouds=.false.
   logical :: do_donner_deep=.false.
   logical :: do_dryadj=.false.
   logical :: do_lin_cld_microphys=.false.
   logical :: do_lsc=.true.
   logical :: do_mca=.true. 
   logical :: do_ras=.false.
   logical :: do_rh_clouds=.false.
   logical :: do_strat=.false.
   logical :: do_uw_conv=.false.
! tracers 
   logical :: do_tracers_in_donner =.false.
   logical :: do_tracers_in_mca = .false.
   logical :: do_tracers_in_ras = .false.
   logical :: do_tracers_in_uw = .false.
! donner specific 
   logical :: do_donner_before_uw = .false.
   logical :: do_donner_mca=.true.
   logical :: do_donner_conservation_checks = .false.
   logical :: do_limit_donner = .false. ! .false. produces previous 
                                        ! behavior (cjg)
   logical :: force_donner_moist_conserv = .false.
! cmt specific
   logical :: cmt_uses_donner = .false.
   logical :: cmt_uses_ras = .false.
   logical :: cmt_uses_uw  = .false.
! others
   logical :: doing_diffusive
   logical :: use_updated_profiles_for_uw = .false.
   logical :: only_one_conv_scheme_per_column = .false.
   logical :: limit_conv_cloud_frac = .false.

! ---> h1g
   real    :: conv_frac_max = 0.99
   logical :: use_updated_profiles_for_clubb = .false.
   logical :: remain_detrain_bug = .false.
! <--- h1g

   logical :: include_donmca_in_cosp = .true.
   logical :: use_tau=.false.
   logical :: do_gust_cv = .false.
   logical :: do_liq_num = .false.
   logical :: do_simple =.false.
   logical :: do_unified_convective_closure = .false.
   logical :: do_limit_uw = .false.     ! .false. produces previous
                                        ! behavior (cjg )
   logical :: using_fms = .true.
   logical :: do_ice_num=.false.
   logical :: detrain_liq_num=.false.
   logical :: detrain_ice_num =.false.
   logical :: do_legacy_strat_cloud = .true.
   character(len=64)  :: cmt_mass_flux_source = 'ras'

   integer :: tau_sg = 0
   integer :: k_sg = 2

   real :: pdepth = 150.e2
   real :: gustmax = 3.                    ! maximum gustiness wind (m/s)
   real :: gustconst = 10./SECONDS_PER_DAY ! constant in kg/m2/sec, default =
                                           ! 1 cm/day = 10 mm/day

namelist /moist_processes_nml/ do_mca, do_lsc, do_ras, do_uw_conv, do_strat,     &
                               do_donner_before_uw, use_updated_profiles_for_uw, &
                               only_one_conv_scheme_per_column, do_diag_clouds,  &
                               limit_conv_cloud_frac, do_dryadj, pdepth,         &
                               include_donmca_in_cosp, &
                               do_unified_convective_closure, tau_sg, k_sg,      &
                               do_lin_cld_microphys, use_tau, do_rh_clouds,      &
                               cmt_mass_flux_source, do_donner_deep, do_cmt,     &
                               do_gust_cv, cmt_mass_flux_source, gustmax,        &
                               gustconst, do_liq_num, force_donner_moist_conserv,&
                               do_donner_conservation_checks, do_donner_mca,     &
                               do_limit_uw, do_limit_donner, using_fms,          &
                               do_bm, do_bmmass, do_bmomp, do_simple, &
                               do_ice_num, do_legacy_strat_cloud, &
                               detrain_liq_num, detrain_ice_num,  &
                               conv_frac_max, use_updated_profiles_for_clubb, remain_detrain_bug !h1g

!-------------------- clock definitions --------------------------------

integer :: convection_clock, largescale_clock, donner_clock, mca_clock, ras_clock, &
           donner_mca_clock, bm_clock, cmt_clock, closure_clock, lscalecond_clock, &
           stratcloud_clock, shallowcu_clock

!-------------------- diagnostics fields -------------------------------
! ---> h1g, dump cell and neso cloud fraction from donner-deep, 2011-08-08
integer :: id_cell_cld_frac,  id_meso_cld_frac, id_donner_humidity_area
! <--- h1g, dump cell and neso cloud fraction from donner-deep, 2011-08-08

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_snow_tot, id_tot_cld_amt, id_conv_freq, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_LWP, id_IWP, id_AWP, id_gust_conv, &

           id_tot_cloud_area,  id_tot_liq_amt,  id_tot_ice_amt,  &
           id_tot_h2o, id_tot_vapor, &
           id_lsc_cloud_area,  id_lsc_liq_amt,  id_lsc_ice_amt,  &
           id_conv_cloud_area, id_conv_liq_amt, id_conv_ice_amt, &
           id_LWP_all_clouds,  id_IWP_all_clouds, id_WP_all_clouds, &

           id_tdt_dadj, id_rh,  id_qs, id_mc, id_mc_donner, id_mc_full, &
           id_mc_donner_half, &
           id_rh_cmip, id_mc_conv_up, id_mc_half, &
           id_conv_cld_base, id_conv_cld_top, &
           id_tdt_deep_donner, id_qdt_deep_donner, &
           id_qadt_deep_donner, id_qldt_deep_donner, &
           id_qidt_deep_donner, &
           id_qndt_deep_donner,  id_qnidt_deep_donner, &
           id_tdt_mca_donner, id_qdt_mca_donner, &
           id_prec_deep_donner, id_prec_mca_donner,&
           id_tdt_uw, id_qdt_uw, &
           id_qadt_uw, id_qldt_uw, id_qidt_uw, id_qndt_uw, id_qnidt_uw, &
           id_prec1_deep_donner, &
           id_snow_deep_donner, id_snow_mca_donner, &
           id_qadt_ls, id_qldt_ls, id_qndt_ls, id_qidt_ls, id_qnidt_ls, &
           id_qadt_conv, id_qldt_conv, id_qndt_conv, id_qidt_conv, &
           id_qnidt_conv, &
           id_qa_ls_col, id_ql_ls_col, id_qn_ls_col, id_qi_ls_col, &
           id_qni_ls_col, &
           id_qa_conv_col, id_ql_conv_col, id_qn_conv_col,  &
           id_qni_conv_col, id_qi_conv_col, &
           id_bmflag, id_klzbs, id_invtaubmt, id_invtaubmq, &
           id_massflux, id_entrop_ls, &
           id_cape, id_cin, id_tref, id_qref, &
           id_q_conv_col, id_q_ls_col, id_t_conv_col, id_t_ls_col, &
           id_enth_moist_col, id_wat_moist_col, &
           id_enth_ls_col, id_wat_ls_col, &
           id_enth_conv_col, id_wat_conv_col, &
           id_enth_donner_col, id_wat_donner_col, &
           id_enth_donner_col2,  &
           id_enth_donner_col3,  &
           id_enth_donner_col4,  &
           id_enth_donner_col5,  &
           id_enth_donner_col6,  &
           id_enth_donner_col7,  &
           id_enth_mca_donner_col, id_wat_mca_donner_col, &
           id_enth_uw_col, id_wat_uw_col, &
           id_scale_donner, id_scale_uw, &
           id_ras_precip, id_ras_freq, id_don_precip, id_don_freq, &
           id_lsc_precip, id_lsc_freq, id_uw_precip, id_uw_snow, &
           id_uw_freq, &
           id_prod_no, id_m_cdet_donner, id_m_cellup, &
           id_conv_rain3d, id_conv_snow3d,   &
           id_lscale_rain3d, id_lscale_snow3d, id_lscale_precip3d
 
integer :: id_qvout, id_qaout, id_qlout, id_qiout
integer :: id_qnout, id_qniout

integer :: id_vaporint, id_condensint, id_precipint, id_diffint
integer :: id_vertmotion
integer :: id_max_enthalpy_imbal_don, id_max_water_imbal_don
integer :: id_max_enthalpy_imbal, id_max_water_imbal
integer :: id_enthint, id_lprcp, id_lcondensint, id_enthdiffint
integer :: id_wetdep_om, id_wetdep_SOA, id_wetdep_bc, &
           id_wetdep_so4, id_wetdep_so2, id_wetdep_DMS, &
           id_wetdep_NH4NO3, id_wetdep_salt, id_wetdep_dust
integer :: id_f_snow_berg, id_f_snow_berg_cond, id_f_snow_berg_wtd

integer, dimension(:), allocatable :: id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col, &
                                      id_tracerdt_mcadon, &
                                      id_tracerdt_mcadon_col, &
                                      id_wetdep, &
                                      id_wet_deposition
real :: missing_value = -999.

!-------------------- individual scheme tracers ------------------------
   logical, dimension(:), allocatable :: tracers_in_donner, tracers_in_uw, &
                                         tracers_in_mca, tracers_in_ras
   integer :: num_donner_tracers=0
   integer :: num_mca_tracers=0
   integer :: num_ras_tracers=0
   integer :: num_uw_tracers=0
   integer :: num_tracers=0

   integer :: nbcphobic =0
   integer :: nbcphilic =0
   integer :: nomphobic =0
   integer :: nomphilic =0
   integer :: nsalt1 =0
   integer :: nsalt2 =0
   integer :: nsalt3 =0
   integer :: nsalt4 =0
   integer :: nsalt5 =0
   integer :: ndust1    =0
   integer :: ndust2    =0
   integer :: ndust3    =0
   integer :: ndust4    =0
   integer :: ndust5    =0
   integer :: nDMS      =0
   integer :: nSO2      =0
   integer :: nSO4      =0
   integer :: nSOA      =0
   integer :: nNH4NO3   =0
   integer :: nNH4      =0
   

!------------------- other global variables and parameters -------------
   real, parameter :: epst=200.

   integer :: nsphum, nql, nqi, nqa, nqn   ! tracer indices for stratiform clouds
   integer :: nqni
   integer :: nqr, nqs, nqg                ! additional tracer indices for Lin Micro-Physics
   integer :: ktop                         ! top layer index for Lin Micro-Physics
   logical :: do_cosp, donner_meso_is_largescale
   real    :: strat_precip_in_cosp = 0.
   real    :: donner_precip_in_cosp = 0.
   real    :: uw_precip_in_cosp = 0.
!-->cjg
   integer :: do_clubb
!<--cjg


!------------------ allocatable moist processes variables --------------

   real, allocatable, dimension(:,:)   :: max_enthalpy_imbal, max_water_imbal, &
                                          max_enthalpy_imbal_don, max_water_imbal_don
   real, allocatable, dimension(:,:,:) :: tin, qin, rin, uin, vin, &
                                          ttnd, qtnd, rtnd, utnd, vtnd, ttnd_don, qtnd_don, &
                                          delta_temp, delta_vapor, delta_q, &
                                          donner_humidity_area, donner_humidity_factor
   real, allocatable, dimension(:,:,:) :: delta_qni, delta_qn
   real, allocatable, dimension(:,:,:) :: nllin, nilin
   real, allocatable, dimension(:,:,:) :: tin_orig, qin_orig, tdt_init, qdt_init
   real, allocatable, dimension(:,:,:) :: qtnd_wet,  &         ! specific humidity tendency (kg/kg/s)
                                          cloud_wet, &         ! cloud liquid+ice (kg/kg)
                                          cloud_frac           ! cloud area fraction
   real, allocatable, dimension(:,:,:) :: liquid_precip, frozen_precip
   real, allocatable, dimension(:,:,:) :: frz_meso, liq_meso, frz_cell
   real, allocatable, dimension(:,:,:) :: liq_cell, mca_frz, mca_liq
   real, allocatable, dimension(:,:,:) :: frz_mesoh, liq_mesoh, frz_cellh, &
                                          liq_precflx, ice_precflx, &
                                          liq_cellh, mca_frzh, mca_liqh,&
                                          ice_precflxh, liq_precflxh
   real, allocatable, dimension(:,:) ::   sumneg
   real, allocatable, dimension(:,:,:) :: ttnd_conv, qtnd_conv
   real, allocatable, dimension(:,:,:) :: qsat, det0, det_cmt       
   real, allocatable, dimension(:,:,:) :: mc_full, mc_donner, m_cdet_donner, massflux, mc_donner_up, &
                                          mc_half, mc_donner_half
   real, allocatable, dimension(:,:,:) :: RH, wetdeptnd, q_ref, t_ref
   real, allocatable, dimension(:,:,:) :: cf, cmf
   real, allocatable, dimension(:,:,:,:) :: tracer,tracer_orig, rdt_init, &
                                            qtr, q_tnd, donner_tracer

   real, allocatable, dimension(:,:)   :: prec_intgl  

! ---> h1g, save cloud condensate tendency due to convection (20120817) 
   real, allocatable, dimension(:,:,:) :: qldt_conv, qidt_conv, qadt_conv, qndt_conv, qnidt_conv
! <--- h1g
!-----------------------------------------------------------------------

                             contains

!#######################################################################
! used to allocate variables used throughout moist_processes
!--> cjg: code modification to allow diagnostic tracers in physics_up (20120508) 
!         lx is the number of prognostic tracers
!         mx is the total number of tracers (prognostic+diagnostic)

!subroutine moist_alloc_init (ix, jx, kx, lx)
!   integer, intent(in) :: ix,jx,kx,lx

subroutine moist_alloc_init (ix, jx, kx, lx, mx)
   integer, intent(in) :: ix,jx,kx,lx,mx
!<--cjg

   if (moist_allocated) return

   allocate( tin       (ix,jx,kx))                          !; tin                    = 0.0
   allocate( qin       (ix,jx,kx))                          !; qin                    = 0.0
   allocate( rin       (ix,jx,kx))                          !; rin                    = 0.0
   allocate( uin       (ix,jx,kx))                          !; uin                    = 0.0
   allocate( vin       (ix,jx,kx))                          !; vin                    = 0.0
   allocate( tin_orig  (ix,jx,kx))                          !; tin_orig               = 0.0
   allocate( qin_orig  (ix,jx,kx))                          !; qin_orig               = 0.0
   allocate( t_ref     (ix,jx,kx))                          ; t_ref                  = 0.0
   allocate( q_ref     (ix,jx,kx))                          ; q_ref                  = 0.0
   allocate( ttnd      (ix,jx,kx))                          ; ttnd                   = 0.0
   allocate( qtnd      (ix,jx,kx))                          ; qtnd                   = 0.0
   allocate( rtnd      (ix,jx,kx))                          ; rtnd                   = 0.0
   allocate( utnd      (ix,jx,kx))                          ; utnd                   = 0.0
   allocate( vtnd      (ix,jx,kx))                          ; vtnd                   = 0.0
   allocate( ttnd_don  (ix,jx,kx))                          ; ttnd_don               = 0.0
   allocate( qtnd_don  (ix,jx,kx))                          ; qtnd_don               = 0.0
   allocate( ttnd_conv (ix,jx,kx))                          ; ttnd_conv              = 0.0
   allocate( qtnd_conv (ix,jx,kx))                          ; qtnd_conv              = 0.0
   allocate( qtnd_wet  (ix,jx,kx))                          ; qtnd_wet               = 0.0
   allocate( tdt_init  (ix,jx,kx))                          ; tdt_init               = 0.0
   allocate( qdt_init  (ix,jx,kx))                          ; qdt_init               = 0.0
   allocate( cf        (ix,jx,kx))                          ; cf                     = 0.0
   allocate( cmf       (ix,jx,kx))                          ; cmf                    = 0.0
   allocate( delta_temp(ix,jx,kx))                          ; delta_temp             = 0.0
   allocate( delta_q   (ix,jx,kx))                          ; delta_q                = 0.0
   allocate( delta_vapor(ix,jx,kx))                         ; delta_vapor            = 0.0
   allocate( donner_humidity_area(ix,jx,kx))                ; donner_humidity_area   = 0.0
   allocate( donner_humidity_factor(ix,jx,kx))              ; donner_humidity_factor = 0.0
   allocate( cloud_wet  (ix,jx,kx))                         ; cloud_wet              = 0.0
   allocate( cloud_frac (ix,jx,kx))                         ; cloud_frac             = 0.0
   allocate( liquid_precip(ix,jx,kx))                       ; liquid_precip          = 0.0
   allocate( frozen_precip(ix,jx,kx))                       ; frozen_precip          = 0.0
   allocate( ice_precflx (ix,jx,kx))                        ; ice_precflx            = 0.0
   allocate( liq_precflx (ix,jx,kx))                        ; liq_precflx            = 0.0
   allocate( frz_meso  (ix,jx,kx))                          ; frz_meso               = 0.0
   allocate( liq_meso  (ix,jx,kx))                          ; liq_meso               = 0.0
   allocate( frz_cell  (ix,jx,kx))                          ; frz_cell               = 0.0
   allocate( liq_cell  (ix,jx,kx))                          ; liq_cell               = 0.0
   allocate( mca_frz   (ix,jx,kx))                          ; mca_frz                = 0.0
   allocate( mca_liq   (ix,jx,kx))                          ; mca_liq                = 0.0
   allocate( frz_mesoh (ix,jx,kx+1))                        ; frz_mesoh              = 0.0
   allocate( liq_mesoh (ix,jx,kx+1))                        ; liq_mesoh              = 0.0
   allocate( frz_cellh (ix,jx,kx+1))                        ; frz_cellh              = 0.0
   allocate( sumneg    (ix,jx))                        ; sumneg                 = 0.0
   allocate( liq_cellh (ix,jx,kx+1))                        ; liq_cellh              = 0.0
   allocate( mca_liqh  (ix,jx,kx+1))                        ; mca_liqh               = 0.0
   allocate( mca_frzh  (ix,jx,kx+1))                        ; mca_frzh               = 0.0
   allocate( ice_precflxh(ix,jx,kx+1))                      ; ice_precflxh           = 0.0
   allocate( liq_precflxh(ix,jx,kx+1))                      ; liq_precflxh           = 0.0
   allocate( qsat      (ix,jx,kx))                          ; qsat                   = 0.0
   allocate( det0      (ix,jx,kx))                          ; det0                   = 0.0
   allocate( det_cmt   (ix,jx,kx))                          ; det_cmt                = 0.0
   allocate( mc_full   (ix,jx,kx))                          ; mc_full                = 0.0
   allocate( mc_donner (ix,jx,kx))                          ; mc_donner              = 0.0
   allocate( mc_donner_up (ix,jx,kx))                       ; mc_donner_up           = 0.0
   allocate( mc_half      (ix,jx,kx+1))                     ; mc_half                = 0.0
   allocate( mc_donner_half (ix,jx,kx+1))                   ; mc_donner_half         = 0.0
   allocate( m_cdet_donner(ix,jx,kx))                       ; m_cdet_donner          = 0.0
   allocate( massflux  (ix,jx,kx))                          ; massflux               = 0.0
   allocate( RH        (ix,jx,kx))                          ; RH                     = 0.0
! pmass defined in moist_processes_utils
   allocate( pmass     (ix,jx,kx))                          ; pmass                  = 0.0
   allocate( wetdeptnd (ix,jx,kx))                          ; wetdeptnd              = 0.0
!--> cjg: code modification to allow diagnostic tracers in physics_up (20120508) 
!   allocate(tracer     (ix,jx,kx,lx))                       ; tracer                 = 0.0
!   allocate(tracer_orig(ix,jx,kx,lx))                       ; tracer_orig            = 0.0
   allocate(tracer     (ix,jx,kx,mx))                       ; tracer                 = 0.0
   allocate(tracer_orig(ix,jx,kx,mx))                       ; tracer_orig            = 0.0
!<--cjg
   allocate(q_tnd      (ix,jx,kx,lx))                       ; q_tnd                  = 0.0
   allocate(rdt_init   (ix,jx,kx,lx))                       ; rdt_init               = 0.0
   allocate(qtr          (ix,jx,kx,num_donner_tracers))     ; qtr                    = 0.0
   allocate(donner_tracer(ix,jx,kx,num_donner_tracers))     ; donner_tracer          = 0.0
   allocate(delta_qn   (ix,jx,kx))                          ; delta_qn               = 0.0
   allocate(delta_qni  (ix,jx,kx))                          ; delta_qni              = 0.0
   allocate(nllin      (ix,jx,kx))                          ; nllin                  = 0.0
   allocate(nilin      (ix,jx,kx))                          ; nilin                  = 0.0

! ---> h1g, allocate cloud condensate tendency due to convection (20120817) 
   allocate( qldt_conv (ix,jx,kx))
   allocate( qidt_conv (ix,jx,kx))
   allocate( qadt_conv (ix,jx,kx))
   if( do_liq_num ) allocate( qndt_conv (ix,jx,kx))
   if( do_ice_num ) allocate( qnidt_conv (ix,jx,kx))
! <--- h1g

   moist_allocated = .true.
  
end subroutine moist_alloc_init


!#######################################################################
! used to deallocate variables used throughout moist_processes
subroutine moist_alloc_end

   if (moist_allocated .eqv. .false. ) return
   deallocate( tin       )
   deallocate( qin       )
   deallocate( rin       )
   deallocate( uin       )
   deallocate( vin       )
   deallocate( tin_orig  )
   deallocate( qin_orig  )
   deallocate( t_ref     )
   deallocate( q_ref     )
   deallocate( ttnd      )
   deallocate( qtnd      )
   deallocate( rtnd      )
   deallocate( utnd      )
   deallocate( vtnd      )
   deallocate( ttnd_don  )
   deallocate( qtnd_don  )
   deallocate( ttnd_conv )
   deallocate( qtnd_conv )
   deallocate( qtnd_wet  )
   deallocate( tdt_init  )
   deallocate( qdt_init  )
   deallocate( cf        )
   deallocate( cmf       )
   deallocate( delta_temp)
   deallocate( delta_q   )
   deallocate( delta_vapor )
   deallocate( donner_humidity_area)
   deallocate( donner_humidity_factor)
   deallocate( cloud_wet  )
   deallocate( cloud_frac )
   deallocate( liquid_precip)
   deallocate( frozen_precip)
   deallocate( ice_precflx)
   deallocate( liq_precflx)
   deallocate( frz_meso  )
   deallocate( liq_meso  )
   deallocate( frz_cell  )
   deallocate( liq_cell  )
   deallocate( mca_frz   )
   deallocate( mca_liq   )
   deallocate( frz_mesoh )
   deallocate( liq_mesoh )
   deallocate( frz_cellh )
   deallocate( sumneg    )
   deallocate( liq_cellh )
   deallocate( mca_frzh  )
   deallocate( mca_liqh  )
   deallocate( ice_precflxh)
   deallocate( liq_precflxh)
   deallocate( qsat      )
   deallocate( det0      )
   deallocate( det_cmt   )
   deallocate( mc_full   )
   deallocate( mc_donner )
   deallocate( mc_donner_up )
   deallocate( mc_half      )
   deallocate( mc_donner_half      )
   deallocate( m_cdet_donner)
   deallocate( massflux  )
   deallocate( RH        )
   deallocate( pmass     )
   deallocate( wetdeptnd )
   deallocate(tracer     )
   deallocate(tracer_orig)
   deallocate(q_tnd      )
   deallocate(rdt_init   )
   deallocate(qtr        )
   deallocate(donner_tracer)
   deallocate(delta_qn   )
   deallocate(delta_qni  )
   deallocate(nllin      )
   deallocate(nilin      )

! ---> h1g, deallocate cloud condensate tendency due to convection (20120817) 
   deallocate( qldt_conv )
   deallocate( qidt_conv )
   deallocate( qadt_conv )
   if( do_liq_num ) deallocate( qndt_conv )
   if( do_ice_num ) deallocate( qnidt_conv )
! <--- h1g

   moist_allocated = .false.

end subroutine moist_alloc_end

!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            radturbten, cush, cbmf,                    &
                            pblht, ustar, bstar, qstar,                &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
                            convect, lprec, fprec, fl_lsrain,          &
                            fl_lssnow, fl_ccrain, fl_ccsnow, &
                            fl_donmca_rain, fl_donmca_snow, gust_cv,  &
                            area, lon, lat,                           &
                            mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land       fraction of surface covered by land
!                      [real, dimension(nlon,nlat)]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega      omega (vertical velocity) at full levels
!                    in pascals per second
!                      [real, dimension(nlon,nlat,nlev)]
!
!         diff_t     vertical diffusion coefficient for temperature
!                    and tracer (m*m/sec) on half levels
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels,
!                    at the current time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
! 
!         tm, qm     temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the previous time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rm         tracer fields at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         um, vm     zonal and meridional wind [m/s] at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev)]
!
!         area       grid box area (in m2)
!                      [real, dimension(nlon,nlat)]
!
!         lon        longitude in radians           ! h1g
!                      [real, dimension(nlon,nlat)] ! h1g
!
!         lat        latitude in radians
!                      [real, dimension(nlon,nlat)]
!  
! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]
! 
!   out:  convect    is moist convection occurring in this grid box?
!                      [logical, dimension(nlon,nlat)]
!
!         lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
! 
!         gust_cv    gustiness from convection  in m/s
!                      [real, dimension(nlon,nlat)]
!
!       optional
!  -----------------
! 
!    in:  mask       mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!         kbot       index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!
!-----------------------------------------------------------------------
   integer,         intent(in)           :: is,ie,js,je
   type(time_type), intent(in)           :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:)     :: land, pblht, ustar, bstar, qstar
   real, intent(inout), dimension(:,:)   :: cush, cbmf
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull, zhalf, zfull, omega, &
                                            diff_t, t, q, u, v, tm, qm, um, vm
   real, dimension(:,:,:), intent(in)    :: radturbten
   real, intent(in), dimension(:,:,:,:) :: r, rm                      ! cjg: inout
   real, intent(inout),dimension(:,:,:)  :: tdt, qdt, udt, vdt
   real, intent(inout),dimension(:,:,:,:):: rdt
logical, intent(out), dimension(:,:)     :: convect
   real, intent(out), dimension(:,:)     :: lprec, fprec, gust_cv
   real, intent(out), dimension(:,:,:)   :: fl_lsrain, fl_lssnow, &
                                            fl_ccrain, fl_ccsnow, &
                                            fl_donmca_rain, fl_donmca_snow
   real, intent(out), dimension(:,:,:)   :: diff_cu_mo
   real, intent(in) , dimension(:,:)     :: area
   real, intent(in) , dimension(:,:)     :: lon
   real, intent(in) , dimension(:,:)     :: lat


   real, intent(in) , dimension(:,:,:), optional :: mask
   integer, intent(in), dimension(:,:), optional :: kbot

!-----------------------------------------------------------------------
   integer :: secs, days
   integer :: n, nn, i, j, k, ix, jx, kx, nt, tr
   integer :: m, mm
   logical :: used, avgbl
   real    :: dtinv

   real, dimension(size(t,1),size(t,2)) :: cape, cin
   real, dimension(size(t,1),size(t,2)) :: precip, total_precip, lheat_precip, &
                                           precip_returned, precip_adjustment, &
                                           vert_motion
   real, dimension(size(t,1),size(t,2)) :: rain, snow, &
                                           rain_don, snow_don, &
                                           rain_ras, snow_ras, &
                                           rain_donmca, snow_donmca
   real, dimension(size(t,1),size(t,2)) :: bmflag, klzbs, invtaubmt, invtaubmq
   real, dimension(size(t,1),size(t,2)) :: scale
   real, dimension(size(t,1),size(t,2)) :: freq_count
   real, dimension(size(t,1),size(t,2)) :: enthint, lcondensint, enthdiffint,  &
                                           vaporint, condensint, precipint, diffint

   real, dimension(size(t,1),size(t,2),size(phalf,3)) :: rain3d, snow3d
   real, dimension(size(t,1),size(t,2),size(phalf,3)) :: snowclr3d
   real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc, m_cellup, mc_cmt
   real, dimension(size(t,1),size(t,2),size(pfull,3)) :: f_snow_berg


!     sfc_sh_flux      sensible heat flux across the surface
!                      [ watts / m**2 ]
!     sfc_vapor_flux   water vapor flux across the surface
!                      [ kg(h2o) / (m**2 sec) ]
!     tr_flux          tracer fux across the surface
!                      [ kg(tracer) / (m**2 sec) ]
   real, dimension(size(t,1),size(t,2)) :: sfc_sh_flux, sfc_vapor_flux
   real, dimension(size(t,1),size(t,2),num_donner_tracers) :: tr_flux  
   real, dimension(size(t,1),size(t,2),num_donner_tracers) :: &
                                                          donner_wetdep
   real, dimension(size(t,1),size(t,2),num_uw_tracers) :: &
                                                          uw_wetdep
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) :: total_wetdep
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) ::  &
                                                       total_wetdep_uw
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) ::   &
                                                     total_wetdep_donner
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) :: ls_wetdep
   real, dimension(size(t,1),size(t,2),size(t,3) ) :: total_conv_cloud,&
                           conv_cld_frac, tot_conv_liq, tot_conv_ice

!chemistry start
   real, parameter :: boltz = 1.38044e-16
   integer, dimension(size(rdt,1),size(rdt,2)) :: cldtop, cldbot
   real, dimension(size(rdt,1),size(rdt,2),size(rdt,3)) :: prod_no
   real, dimension(size(rdt,1),size(rdt,2),size(rdt,3),size(rdt,4)) :: wet_data
!chemistry end

   real, dimension(size(t,1),size(t,2))           ::  adjust_frac      
   real, dimension(size(t,1),size(t,2),size(t,3)) ::  ttnd_adjustment
   real, dimension(size(t,1),size(t,2),size(t,3)) ::  available_cf_for_uw

   logical, dimension(size(t,1),size(t,2)) :: conv_calc_completed
   logical, dimension(size(t,1),size(t,2)) :: coldT

!temporary variables
   real :: temp
   logical, dimension(size(t,1),size(t,2)) :: ltemp
   real, dimension(size(t,1),size(t,2)) :: temp_2d
   real, dimension(size(t,1),size(t,2)) :: tca2
   real, dimension(size(t,1),size(t,2),size(t,3)) :: total_cloud_area
   real, dimension(size(t,1),size(t,2),size(t,3)) :: temp_3d1, temp_3d2, temp_3d3

! ---> h1g, 2010-08-23
   real                                 ::       current_total_sec
   integer                              ::       current_sec, current_days

!  consider the donner-deep mass flux impacts on clubb
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: conv_frac_clubb
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: convective_humidity_ratio_clubb
   real                                                        :: qrf, env_fraction, env_qv
! <--- h1g, 2010-08-23


! ---> h1g, 2012-10-05
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: qcvar_clubb
! <--- h1g, 2012-10-05      
!-------- input array size and position in global storage --------------
   write(6,*) 'start moist processes 0'

      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(rdt,4)

       
!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('moist_processes_mod',  &
                 'moist_processes_init has not been called.', FATAL)
      endif

      conv_calc_completed = .false.
      available_cf_for_uw = 1.0

!--------------------------------------------------------------------
!    define the inverse of the time step.
!--------------------------------------------------------------------
      dtinv = 1.0/dt

!--------------------------------------------------------------------
!    initialize the arrays which will be used in this subroutine.
!--------------------------------------------------------------------
      write(6,*) 'start moist processes 0.5'

      rain_don     = 0.0
      snow_don     = 0.0
      rain_donmca  = 0.0
      snow_donmca  = 0.0
      write(6,*) 'start moist processes 0.6'

      lprec        = 0.0  
      fprec        = 0.0
      write(6,*) 'start moist processes 0.7'

      write(6,*) size(fl_lsrain,1),size(fl_lsrain,2),size(fl_lsrain,3)
      fl_lsrain(:,:,:) = 0.
      fl_lssnow(:,:,:) = 0.
      fl_ccrain(:,:,:) = 0.
      write(6,*) 'start moist processes 0.75'

      fl_ccsnow(:,:,:) = 0.
      fl_donmca_rain(:,:,:) = 0.
      fl_donmca_snow(:,:,:) = 0.
      convect      = .false.
      gust_cv      = 0.0
      precip       = 0.0 
      rain3d       = 0.0
      snow3d       = 0.0
      write(6,*) 'start moist processes 1', size(rdt,1), size(rdt,2), size(rdt,3), size(rdt,4)
!---------------------------------------------------------------------
!    initialize local arrays which will hold sums.
!---------------------------------------------------------------------
      rdt_init(is:ie,js:je,:,:)  = rdt
      tdt_init(is:ie,js:je,:)  = tdt
      qdt_init(is:ie,js:je,:)  = qdt
!      ttnd_conv(is:ie,js:je,:) = 0.
!      qtnd_conv(is:ie,js:je,:) = 0.
!      qtnd(is:ie,js:je,:)      = 0.
!      q_tnd(is:ie,js:je,:,:)     = 0.

!---------------------------------------------------------------------
!    define input fields to be used, either the tau time level fields,
!    or the tau - 1 time level values updated with the time tendencies
!    thus far calculated on the current step. control is through nml
!    variable use_tau.
!---------------------------------------------------------------------
      write(6,*) 'start moist processes 2'

      if (use_tau) then
        tin(is:ie,js:je,:) = t
        qin(is:ie,js:je,:) = q
        uin(is:ie,js:je,:) = u
        vin(is:ie,js:je,:) = v
        do tr=1,size(r,4)
          tracer(is:ie,js:je,:,tr) = r(:,:,:,tr)
        end do  
      else
        tin(is:ie,js:je,:) = tm + tdt*dt
        qin(is:ie,js:je,:) = qm + qdt*dt
        uin(is:ie,js:je,:) = um + udt*dt
        vin(is:ie,js:je,:) = vm + vdt*dt
        do tr=1,size(rdt,4)
          tracer(is:ie,js:je,:,tr) = rm(:,:,:,tr) + rdt(:,:,:,tr)*dt
        end do  
        do tr=size(rdt,4) +1, size(r,4)
          tracer(is:ie,js:je,:,tr) = r(:,:,:,tr)
        end do  
      endif

!--------------------------------------------------------------------
!    if using eta vertical coordinate, define the appropriate values 
!    for any points located below the ground. values of 0.0 are given
!    to u, v and q, and a temperature value of EPST (=200. K) is given 
!    to sub-surface  points.
!--------------------------------------------------------------------
      write(6,*) 'start moist processes 3'

      if (present(mask) .and. present(kbot))  then
        tin(is:ie,js:je,:) = mask*tin(is:ie,js:je,:) + (1.0 - mask)*EPST 
        qin(is:ie,js:je,:) = mask*qin(is:ie,js:je,:)
        uin(is:ie,js:je,:) = mask*uin(is:ie,js:je,:)
        vin(is:ie,js:je,:) = mask*vin(is:ie,js:je,:)
        do tr=1,size(r,4)
          tracer(is:ie,js:je,:,tr) = mask(:,:,:)*tracer(is:ie,js:je,:,tr)
        end do  
      endif
   
!----------------------------------------------------------------------
!    compute the mass in each model layer.
!----------------------------------------------------------------------
      write(6,*) 'start moist processes 4'

      do k=1,kx
        pmass(is:ie,js:je,k) = (phalf(:,:,k+1) - phalf(:,:,k))/GRAV
      end do

!----------------------------------------------------------------------
!    output any requested convectively-transported tracer fields 
!    and / or their column sums before convective transport.
!----------------------------------------------------------------------
      write(6,*) 'start moist processes 5'

      do n=1,num_tracers
        used = send_data (id_conv_tracer(n), tracer(is:ie,js:je,:,n), Time, &
                          is, js, 1, rmask=mask)
        if (id_conv_tracer_col(n) > 0)  &
          call column_diag(id_conv_tracer_col(n), is, js, Time, &
                           tracer(is:ie,js:je,:,n), 1.0) 
      end do

!----------------------------------------------------------------------
!    compute the mean temperature in the lower atmosphere (the lowest
!    pdepth Pa), to be used to determine whether rain or snow reaches
!    the surface. define a logical variable coldT indicating whether
!    snow or rain falls in the column.
!    ????    SHOULD TIN BE USED RATHER THAN t ??
!----------------------------------------------------------------------
      write(6,*) 'start moist processes 6'

      call tempavg (pdepth, phalf, t, snow, mask)
      coldT = .false.
      where (snow(:,:) <= TFREEZE)
        coldT(:,:) = .true.
      endwhere
      
!---------------------------------------------------------------------
!    begin the clock timing the dry and moist convection parameter-
!    izations.
!---------------------------------------------------------------------
      call mpp_clock_begin (convection_clock)

! subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
!                             phalf, pfull, zhalf, zfull, omega, diff_t, &
!                             radturbten, cush, cbmf,                    &
!                             pblht, ustar, bstar, qstar,                &
!                             t, q, r, u, v, tm, qm, rm, um, vm,         &
!                             tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
!                             convect, lprec, fprec, fl_lsrain,          &
!                             fl_lssnow, fl_ccrain, fl_ccsnow, &
!                             fl_donmca_rain, fl_donmca_snow, gust_cv,  &
!                             area, lon, lat, lsc_cloud_area, lsc_liquid,     &
!                             lsc_ice, lsc_droplet_number, &
!                             lsc_ice_number, lsc_snow, lsc_rain,  &
!                             lsc_snow_size, lsc_rain_size     , &
! ! ---> h1g
!                             dcond_ls_liquid,     dcond_ls_ice,         &
!                             Ndrop_act_CLUBB,     Icedrop_act_CLUBB,    &
!                             ndust, rbar_dust,                          &
!                             diff_t_clubb,                              &
!                             tdt_shf,                                   &
!                             qdt_lhf,                                   &
! ! <--- h1g
!                             Aerosol, mask, kbot, &
!                             shallow_cloud_area, shallow_liquid,  &
!                             shallow_ice, shallow_droplet_number, &
!                             shallow_ice_number, &
!                             cell_cld_frac, cell_liq_amt, cell_liq_size, &
!                             cell_ice_amt, cell_ice_size, &
!                             cell_droplet_number, &
!                             meso_cld_frac, meso_liq_amt, meso_liq_size, &
!                             meso_ice_amt, meso_ice_size,  &
!                             meso_droplet_number, nsum_out, &
!                             hydrostatic, phys_hydrostatic)

      write(6,*) ' about to do idm conv and lscale cond'
    call idealized_convection_and_lscale_cond( phalf, pfull, zhalf, zfull, t, r, u, v, tdt, udt, vdt, rdt, Time, dt, mask, kbot)    
    write(6,*) ' Done idm conv and lscale cond'

!---------------------------------------------------------------------
!    end the timing of the convection code section.
!---------------------------------------------------------------------
   call mpp_clock_end (convection_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!              LARGE-SCALE CONDENSATION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!---------------------------------------------------------------------
!    begin the timing of the large-scale condensation code section.
!---------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!         A. NON-PROGNOSTIC CONDENSATION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!    if a non-prognostic cloud scheme is active, then call lscale_cond 
!    to calculate the temperature and specific humidity tendencies 
!    related to the latent heat release associated with the large-scale 
!    supersaturation.
!-----------------------------------------------------------------------
    call mpp_clock_begin (largescale_clock)


!---------------------------------------------------------------------
!    end the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
    call mpp_clock_end (largescale_clock)
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  GENERAL MOISTURE DIAGNOSTICS 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!--------------------------------------------------------------------
!    output diagnostics obtained from the combination of convective and
!    large-scale parameterizations.  
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    output diagnostics obtained from the combination of convective and
!    large-scale parameterizations.  
!--------------------------------------------------------------------


!-----------------------------------------------------------------------
end subroutine moist_processes


!#####################################################################
 
! subroutine moist_processes_time_vary (dt)

! real, intent(in) :: dt


!       if (do_donner_deep) then
!         call donner_deep_time_vary (dt)
!       endif
!       if (do_strat .and. .not. do_legacy_strat_cloud) then
!         call strat_cloud_time_vary (dt, limit_conv_cloud_frac)
!       endif

! end subroutine moist_processes_time_vary


!#####################################################################

subroutine moist_processes_endts (is, js)
 
integer, intent(in) :: is,js

      ! if (do_donner_deep) then
      !   call donner_deep_endts
      ! endif 


      call sum_diag_integral_field ('prec', prec_intgl)
      prec_intgl = 0.0


end subroutine moist_processes_endts



!###################################################################

!#######################################################################
!---> h1g
!subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
subroutine moist_processes_init ( id, jd, kd, lonb, latb, lon, lat, phalf, pref, &
                                  axes, Time, doing_donner)
!<--- h1g

!-----------------------------------------------------------------------
integer,              intent(in)  :: id, jd, kd, axes(4)
real, dimension(:,:), intent(in)  :: lonb, latb
real,dimension(:,:),  intent(in)  :: lon,  lat    ! h1g
real,dimension(:,:,:),intent(in)  :: phalf        ! h1g
real, dimension(:),   intent(in)  :: pref
type(time_type),      intent(in)  :: Time
 logical,              intent(out) :: doing_donner
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer :: unit,io,ierr, n, logunit
character(len=80)  :: scheme
integer            :: secs, days
integer            :: k
!-----------------------------------------------------------------------

       if ( module_is_initialized ) return



       if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
         read (input_nml_file, nml=moist_processes_nml, iostat=io)
         ierr = check_nml_error(io,'moist_processes_nml')
#else

         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'moist_processes_nml')
         enddo
  10     call close_file (unit)
#endif

!--------- write version and namelist to standard log ------------

      call write_version_number(version, tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=moist_processes_nml )

       endif

  
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')
   allocate (prec_intgl(id,jd))


      module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine moist_processes_init

!#######################################################################
!---> h1g
!subroutine moist_processes_end
subroutine moist_processes_end( clubb_term_clock )
integer, intent (out), optional :: clubb_term_clock
!<--- h1g

      if( .not.module_is_initialized ) return


!----------------close various schemes-----------------

! ---> h1g, cjg
      if (do_strat) then
        if (do_clubb > 0) then
! ---> h1g, if CLUBB is in moist-process, CLUBB ends here
          if( do_clubb == 2) then
              call mpp_clock_begin ( clubb_term_clock )
              ! call clubb_end
              call mpp_clock_end ( clubb_term_clock )
          endif
! <--- h1g, if CLUBB is in moist-process, CLUBB ends here
          ! call MG_microp_3D_end
        else
          ! call strat_cloud_end
        end if
      end if
! <--- h1g, cjg

      ! call  detr_ice_num_end
      ! if (do_rh_clouds)   call   rh_clouds_end
      ! if (do_diag_clouds) call  diag_cloud_end
      ! if (do_donner_deep) call donner_deep_end
      ! if (do_cmt        ) call cu_mo_trans_end
      ! if (do_ras        ) call         ras_end
      ! if (do_uw_conv    ) call     uw_conv_end
      ! if (do_lin_cld_microphys) call lin_cld_microphys_end

      deallocate (max_water_imbal)
      deallocate (max_enthalpy_imbal)
      if (do_donner_deep .and. do_donner_conservation_checks) then
        deallocate (max_water_imbal_don)
        deallocate (max_enthalpy_imbal_don)
      endif

      module_is_initialized = .false.

!-----------------------------------------------------------------------

end subroutine moist_processes_end


!#######################################################################
! <SUBROUTINE NAME="moist_processes_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine moist_processes_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

  ! if (do_strat)       call strat_cloud_restart(timestamp)
  ! if (do_diag_clouds) call diag_cloud_restart(timestamp)
  ! if (do_donner_deep) call donner_deep_restart(timestamp)

end subroutine moist_processes_restart
! </SUBROUTINE> NAME="moist_processes_restart"


!#######################################################################

subroutine diag_field_init ( axes, Time, num_tracers, num_donner_tracers )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time
  integer, intent(in) :: num_donner_tracers
  integer, intent(in) :: num_tracers

  character(len=32) :: tracer_units, tracer_name
  character(len=128) :: diaglname
  integer, dimension(3) :: half = (/1,2,4/)
  integer   :: n, nn

!------------ initializes diagnostic fields in this module -------------

   if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then
      id_qref = register_diag_field ( mod_name, &
        'qref', axes(1:3), Time, &
        'Adjustment reference specific humidity profile', &
        'kg/kg',  missing_value=missing_value               )

      id_tref = register_diag_field ( mod_name, &
        'tref', axes(1:3), Time, &
        'Adjustment reference temperature profile', &
        'K',  missing_value=missing_value                   )

      id_bmflag = register_diag_field (mod_name, &
         'bmflag', axes(1:2), Time, &
         'Betts-Miller flag', &
         'no units', missing_value=missing_value            )

      id_klzbs  = register_diag_field  (mod_name, &
         'klzbs', axes(1:2), Time, &
         'klzb', &
         'no units', missing_value=missing_value            )

      id_cape = register_diag_field ( mod_name, & 
        'cape', axes(1:2), Time, &
        'Convectively available potential energy',      'J/Kg')

      id_cin = register_diag_field ( mod_name, &
        'cin', axes(1:2), Time, &
        'Convective inhibition',                        'J/Kg')
   endif

   if ( do_bm ) then
      id_invtaubmt  = register_diag_field  (mod_name, &
         'invtaubmt', axes(1:2), Time, &
         'Inverse temperature relaxation time', &
         '1/s', missing_value=missing_value            )

      id_invtaubmq = register_diag_field  (mod_name, &
         'invtaubmq', axes(1:2), Time, &
         'Inverse humidity relaxation time', &
         '1/s', missing_value=missing_value            )
   end if  ! if ( do_bm )

   if (do_bmmass) then
      id_massflux = register_diag_field (mod_name, &
         'massflux', axes(1:3), Time, &
         'Massflux implied by temperature adjustment', &
         'm/s', missing_value=missing_value                 )
   end if  ! if ( do_bmmass )

   id_ras_precip = register_diag_field ( mod_name, &
     'ras_precip', axes(1:2), Time, &
    'Precipitation rate from ras ',       'kg/m2/s' )

   id_ras_freq = register_diag_field ( mod_name, &
     'ras_freq', axes(1:2), Time, &
    'frequency of precip from ras ',       'number' , &
         missing_value = missing_value                       )

   id_don_precip = register_diag_field ( mod_name, &
     'don_precip', axes(1:2), Time, &
    'Precipitation rate from donner ',       'kg/m2/s' )

   id_don_freq = register_diag_field ( mod_name, &
     'don_freq', axes(1:2), Time, &
    'frequency of precip from donner ',       'number', &
         missing_value = missing_value                       )

   id_lsc_precip = register_diag_field ( mod_name, &
     'lsc_precip', axes(1:2), Time, &
    'Precipitation rate from lsc ',       'kg/m2/s' )

   id_lsc_freq = register_diag_field ( mod_name, &
     'lsc_freq', axes(1:2), Time, &
    'frequency of precip from lsc ',       'number' , &
         missing_value = missing_value                       )

   id_uw_precip = register_diag_field ( mod_name, &
     'uw_precip', axes(1:2), Time, &
    'Precipitation rate from uw shallow',       'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_uw_snow = register_diag_field ( mod_name, &
     'uw_snow', axes(1:2), Time, &
    'Snow rate from uw shallow',       'kg/m2/s' , &
     interp_method = "conserve_order1" )

   id_uw_freq = register_diag_field ( mod_name, &
     'uw_freq', axes(1:2), Time, &
    'frequency of precip from uw shallow ',       'number' , &
         missing_value = missing_value                       )

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from convection ',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from convection ',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from convection ',   'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from convection ','W/m2' )
   
   id_enth_conv_col = register_diag_field ( mod_name, &
     'enth_conv_col', axes(1:2), Time, &
     'Column enthalpy tendency from convection','W/m2' )
 
   id_wat_conv_col = register_diag_field ( mod_name, &
     'wat_conv_col', axes(1:2), Time, &
     'Column total water tendency from convection','kg(h2o)/m2/s' )

   id_enth_donner_col2 = register_diag_field ( mod_name, &
     'enth_donner_col2', axes(1:2), Time, &
     'column enthalpy tendency from Donner liq precip','W/m2' )
 
   id_enth_donner_col3 = register_diag_field ( mod_name, &
     'enth_donner_col3', axes(1:2), Time, &
      'Column enthalpy tendency from Donner frzn precip','W/m2' )
 
   id_enth_donner_col4 = register_diag_field ( mod_name, &
      'enth_donner_col4', axes(1:2), Time, &
     'Atmospheric column enthalpy tendency from Donner convection', &
                                                            'W/m2' )
 
   id_enth_donner_col5 = register_diag_field ( mod_name, &
      'enth_donner_col5', axes(1:2), Time, &
      'Column enthalpy tendency due to condensate xfer from Donner &
                                                &to lsc','W/m2' )

  id_enth_donner_col6 = register_diag_field ( mod_name, &
     'enth_donner_col6', axes(1:2), Time, &
      'Column enthalpy tendency from donner moisture  &
                     &conservation  adjustment','W/m2' )
 
   id_enth_donner_col7 = register_diag_field ( mod_name, &
      'enth_donner_col7', axes(1:2), Time, &
      'Precip adjustment needed to balance donner moisture  &
                                           &adjustment','kg(h2o)/m2/s' )

   id_enth_donner_col = register_diag_field ( mod_name, &
     'enth_donner_col', axes(1:2), Time, &
     'Column enthalpy imbalance from Donner convection','W/m2' )

   id_wat_donner_col = register_diag_field ( mod_name, &
     'wat_donner_col', axes(1:2), Time, &
  'Column total water tendency from Donner convection','kg(h2o)/m2/s' )

   id_enth_mca_donner_col = register_diag_field ( mod_name, &
     'enth_mca_donner_col', axes(1:2), Time, &
    'Column enthalpy imbalance from Donner MCA convection','W/m2' )

   id_wat_mca_donner_col = register_diag_field ( mod_name, &
     'wat_mca_donner_col', axes(1:2), Time, &
     'Column total water imbalance from Donner MCA convection', &
                                                'kg(h2o)/m2/s' )

   id_enth_uw_col = register_diag_field ( mod_name, &
     'enth_uw_col', axes(1:2), Time, &
     'Column enthalpy tendency from UW convection','W/m2' )
 
   id_wat_uw_col = register_diag_field ( mod_name, &
     'wat_uw_col', axes(1:2), Time, &
      'Column total water tendency from UW convection','kg(h2o)/m2/s' )

   id_scale_uw = register_diag_field ( mod_name, &
     'scale_uw', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )
          
   id_scale_donner = register_diag_field ( mod_name, &
     'scale_donner', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_tot  = register_diag_field ( mod_name, &
     'snow_tot ', axes(1:2), Time, &
     'Frozen precip rate from all sources',       'kg(h2o)/m2/s', &
      interp_method = "conserve_order1" )

   id_conv_freq = register_diag_field ( mod_name, &
     'conv_freq', axes(1:2), Time, &
    'frequency of convection ',       'number', &
     missing_value = missing_value                       )

   id_gust_conv = register_diag_field ( mod_name, &
     'gust_conv', axes(1:2), Time, &
    'Gustiness resulting from convection ',       'm/s' )

  id_conv_rain3d= register_diag_field ( mod_name, &
     'conv_rain3d', axes(half), Time, &
    'Rain fall rate from convection -3D ',       'kg(h2o)/m2/s' )

   id_conv_snow3d= register_diag_field ( mod_name, &
     'conv_snow3d', axes(half), Time, &
    'Snow fall rate from convection -3D',       'kg(h2o)/m2/s' )

   id_lscale_rain3d= register_diag_field ( mod_name, &
     'lscale_rain3d', axes(half), Time, &
    'Rain fall rate from lscale  -3D ',   'kg(h2o)/m2/s' )

   id_lscale_snow3d= register_diag_field ( mod_name, &
     'lscale_snow3d', axes(half), Time, &
    'Snow fall rate from lscale -3D',       'kg(h2o)/m2/s' )
   
   id_lscale_precip3d= register_diag_field ( mod_name, &
     'lscale_precip3d', axes(1:3), Time, &
     'LS Precip falling out of gridbox',       'kg(h2o)/m2/s' , &
      mask_variant = .true., missing_value = missing_value)

    id_max_enthalpy_imbal    = register_diag_field    &
       (mod_name, 'max_enth_imbal', axes(1:2), Time,  &
        'max enthalpy  imbalance from moist_processes  ', 'W/m2',   &
              missing_value=missing_value)
    id_max_water_imbal    = register_diag_field    &
         (mod_name, 'max_water_imbal', axes(1:2), Time,   &
      'max water  imbalance from moist_processes  ', 'kg(h2o)/m2/s',  &
              missing_value=missing_value)

    id_enth_moist_col = register_diag_field ( mod_name, &
     'enth_moist_col', axes(1:2), Time, &
     'Column enthalpy imbalance from moist processes','W/m2' )
  
    id_wat_moist_col = register_diag_field ( mod_name, &
      'wat_moist_col', axes(1:2), Time, &
      'Column total water imbalance from moist processes','kg/m2/s' )

    if (do_donner_conservation_checks) then
      id_enthint    = register_diag_field    &
            (mod_name, 'enthint_don', axes(1:2), Time,  &
          'atmospheric column enthalpy change from donner', 'W/m2',  &
          missing_value=missing_value)
     id_lcondensint    = register_diag_field    &
         (mod_name, 'lcondensint_don', axes(1:2), Time, &
         'enthalpy transferred by condensate from donner to lscale', &
            'W/m2',  missing_value=missing_value)
     id_lprcp    = register_diag_field    &
             (mod_name, 'lprcpint_don', axes(1:2),   &
              Time, 'enthalpy removed by donner precip', 'W/m2',   &
             missing_value=missing_value)
      id_vertmotion    = register_diag_field    &
             (mod_name, 'vertmotion_don', axes(1:2), Time,  &
           'enthalpy change due to cell and meso motion in donner',  &
             'W/m2', missing_value=missing_value)
     id_enthdiffint    = register_diag_field    &
            (mod_name, 'enthdiffint_don', axes(1:2),   &
             Time, 'enthalpy  imbalance due to donner', 'W/m2',   &
             missing_value=missing_value)
     id_vaporint    = register_diag_field    &
           (mod_name, 'vaporint_don', axes(1:2),   &
            Time, 'column water vapor change', 'kg(h2o)/m2/s',   &
            missing_value=missing_value)
     id_max_enthalpy_imbal_don    = register_diag_field    &
            (mod_name, 'max_enth_imbal_don', axes(1:2),   &
              Time, 'max enthalpy  imbalance from donner', 'W/m**2',  &
              missing_value=missing_value)
     id_max_water_imbal_don    = register_diag_field    &
            (mod_name, 'max_water_imbal_don', axes(1:2),   &
              Time, 'max water imbalance from donner', 'kg(h2o)/m2/s', &
         missing_value=missing_value)
     id_condensint    = register_diag_field    &
           (mod_name, 'condensint_don', axes(1:2), Time,  &
         'column condensate exported from donner to lscale', &
                         'kg(h2o)/m2/s',  missing_value=missing_value )
     id_precipint    = register_diag_field    &
            (mod_name, 'precipint_don', axes(1:2),   &
             Time, 'column precip from donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
     id_diffint    = register_diag_field    &
          (mod_name, 'diffint_don', axes(1:2),   &
            Time, 'water imbalance due to donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
  endif



if (do_strat ) then

   id_qldt_conv = register_diag_field ( mod_name, &
     'qldt_conv', axes(1:3), Time, &
     'Liquid water tendency from convection',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qndt_conv = register_diag_field ( mod_name, &
     'qndt_conv', axes(1:3), Time, &
     'Liquid drop tendency from convection',      '#/kg/s',  &
                         missing_value=missing_value               )

   id_qidt_conv = register_diag_field ( mod_name, &
     'qidt_conv', axes(1:3), Time, &
     'Ice water tendency from convection',         'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_conv = register_diag_field ( mod_name, &
     'qadt_conv', axes(1:3), Time, &
     'Cloud fraction tendency from convection',    '1/sec',    &
                        missing_value=missing_value               )

   id_ql_conv_col = register_diag_field ( mod_name, &
     'ql_conv_col', axes(1:2), Time, &
    'Liquid water path tendency from convection',  'kg/m2/s' )
   
   id_qn_conv_col = register_diag_field ( mod_name, &
     'qn_conv_col', axes(1:2), Time, &
     'Liquid drp tendency from convection',  'kg/m2/s' )
 
   id_qi_conv_col = register_diag_field ( mod_name, &
     'qi_conv_col', axes(1:2), Time, &
    'Ice water path tendency from convection',     'kg/m2/s' )
   
   id_qa_conv_col = register_diag_field ( mod_name, &
     'qa_conv_col', axes(1:2), Time, &
    'Cloud mass tendency from convection',         'kg/m2/s' )
      
   id_qnidt_conv = register_diag_field ( mod_name, &
     'qnidt_conv', axes(1:3), Time, &
     'Ice number tendency from convection',      '#/kg/s',  &
                         missing_value=missing_value               )

   id_qni_conv_col = register_diag_field ( mod_name, &
     'qni_conv_col', axes(1:2), Time, &
     'Ice number tendency from convection',  'kg/m2/s' )

endif

if ( do_lsc ) then

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
       'Temperature tendency from large-scale cond',   'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from large-scale cond','kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from large-scale cond','W/m2' )
   
 endif

   id_conv_cld_base = register_diag_field ( mod_name, &
     'conv_cld_base', axes(1:2), Time, &
     'pressure at convective cloud base',   'Pa', &
                       mask_variant = .true., &
                       missing_value=missing_value               )

   id_conv_cld_top = register_diag_field ( mod_name, &
     'conv_cld_top', axes(1:2), Time, &
     'pressure at convective cloud top',   'Pa', &
                       mask_variant = .true., &
                       missing_value=missing_value               )

if ( do_strat ) then

   id_mc_full = register_diag_field ( mod_name, &
     'mc_full', axes(1:3), Time, &
     'Net Mass Flux from convection',   'kg/m2/s', &
                       missing_value=missing_value               )
   
   id_mc_half = register_diag_field ( mod_name, &
     'mc_half', axes(half), Time, &
     'Net Mass Flux from convection on half levs',   'kg/m2/s', &
                       missing_value=missing_value               )
   
   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
     'Temperature tendency from strat cloud',        'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from strat cloud',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from strat cloud',          'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from strat cloud',          'kg/m2/s' )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from strat cloud',   'kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from strat cloud','W/m2' )
   
   id_qldt_ls = register_diag_field ( mod_name, &
     'qldt_ls', axes(1:3), Time, &
     'Liquid water tendency from strat cloud',       'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qndt_ls = register_diag_field ( mod_name, &
     'qndt_ls', axes(1:3), Time, &
     'Drop number tendency from strat cloud',        '#/kg/s',  &
                         missing_value=missing_value               )
   id_qidt_ls = register_diag_field ( mod_name, &
     'qidt_ls', axes(1:3), Time, &
     'Ice water tendency from strat cloud',          'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qnidt_ls = register_diag_field ( mod_name, &
     'qnidt_ls', axes(1:3), Time, &
     'Ice number tendency from strat cloud',          '#/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_ls = register_diag_field ( mod_name, &
     'qadt_ls', axes(1:3), Time, &
     'Cloud fraction tendency from strat cloud',     '1/sec',    &
                        missing_value=missing_value               )

   id_ql_ls_col = register_diag_field ( mod_name, &
     'ql_ls_col', axes(1:2), Time, &
    'Liquid water path tendency from strat cloud',   'kg/m2/s' )
   
   id_qn_ls_col = register_diag_field ( mod_name, &
     'qn_ls_col', axes(1:2), Time, &
     'Column drop number tendency from strat cloud',  '#/m2/s' )

   id_qni_ls_col = register_diag_field ( mod_name, &
     'qni_ls_col', axes(1:2), Time, &
     'Column ice particle number tendency from strat cloud',  '#/m2/s' )

   id_qi_ls_col = register_diag_field ( mod_name, &
     'qi_ls_col', axes(1:2), Time, &
    'Ice water path tendency from strat cloud',      'kg/m2/s' )
   
   id_qa_ls_col = register_diag_field ( mod_name, &
     'qa_ls_col', axes(1:2), Time, &
    'Cloud mass tendency from strat cloud',          'kg/m2/s' )
      
   id_enth_ls_col = register_diag_field ( mod_name, &
     'enth_ls_col', axes(1:2), Time, &
     'Column enthalpy tendency from strat cloud','W/m2' )
 
   id_wat_ls_col = register_diag_field ( mod_name, &
     'wat_ls_col', axes(1:2), Time, &
     'Column total water tendency from strat cloud','kg/m2/s' )

endif

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s', &
      interp_method = "conserve_order1" )

   id_WVP = register_diag_field ( mod_name, &
     'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

if ( do_strat ) then

   id_LWP = register_diag_field ( mod_name, &
     'LWP', axes(1:2), Time, &
        'Liquid water path',                            'kg/m2'   )

   id_IWP = register_diag_field ( mod_name, &
     'IWP', axes(1:2), Time, &
        'Ice water path',                               'kg/m2'   )

   id_AWP = register_diag_field ( mod_name, &
     'AWP', axes(1:2), Time, &
        'Column integrated cloud mass ',                'kg/m2'   )

    id_tot_cld_amt = register_diag_field    &
              (mod_name, 'cld_amt_2d', axes(1:2), Time, &
                'total cloud amount', 'percent')

    id_tot_cloud_area = register_diag_field ( mod_name, &
      'tot_cloud_area', axes(1:3), Time, &
      'Cloud area -- all clouds', 'percent', missing_value=missing_value )

    id_tot_h2o     = register_diag_field ( mod_name, &
      'tot_h2o', axes(1:3), Time, &
      'total h2o -- all phases', 'kg/kg', missing_value=missing_value)

    id_tot_vapor     = register_diag_field ( mod_name, &
       'tot_vapor', axes(1:3), Time, &
       'total vapor', 'kg/kg', missing_value=missing_value)

    id_tot_liq_amt = register_diag_field ( mod_name, &
      'tot_liq_amt', axes(1:3), Time, &
      'Liquid amount -- all clouds', 'kg/kg', missing_value=missing_value)

    id_tot_ice_amt = register_diag_field ( mod_name, &
      'tot_ice_amt', axes(1:3), Time, &
      'Ice amount -- all clouds', 'kg/kg', missing_value=missing_value )

    id_lsc_cloud_area = register_diag_field ( mod_name, &
      'lsc_cloud_area', axes(1:3), Time, &
      'Large-scale cloud area', 'percent', missing_value=missing_value )

    id_lsc_liq_amt = register_diag_field ( mod_name, &
      'lsc_liq_amt', axes(1:3), Time, &
      'Large-scale cloud liquid amount', 'kg/kg', missing_value=missing_value )

    id_lsc_ice_amt = register_diag_field ( mod_name, &
      'lsc_ice_amt', axes(1:3), Time, &
      'Large-scale cloud ice amount', 'kg/kg', missing_value=missing_value )

    id_conv_cloud_area = register_diag_field ( mod_name, &
      'conv_cloud_area', axes(1:3), Time, &
      'Convective cloud area', 'percent', missing_value=missing_value )

    id_conv_liq_amt = register_diag_field ( mod_name, &
      'conv_liq_amt', axes(1:3), Time, &
      'Convective cloud liquid amount', 'kg/kg', missing_value=missing_value )

    id_conv_ice_amt = register_diag_field ( mod_name, &
      'conv_ice_amt', axes(1:3), Time, &
      'Convective cloud ice amount', 'kg/kg', missing_value=missing_value)
 
    id_WP_all_clouds = register_diag_field ( mod_name, &
      'WP_all_clouds', axes(1:2), Time, &
      'Total  water path -- all clouds',              'kg/m2'   )

    id_LWP_all_clouds = register_diag_field ( mod_name, &
      'LWP_all_clouds', axes(1:2), Time, &
      'Liquid water path -- all clouds',              'kg/m2'   )

    id_IWP_all_clouds = register_diag_field ( mod_name, &
      'IWP_all_clouds', axes(1:2), Time, &
      'Ice water path -- all clouds',                 'kg/m2'   )

endif

   id_tdt_dadj = register_diag_field ( mod_name, &
     'tdt_dadj', axes(1:3), Time, &
   'Temperature tendency from dry conv adj',       'deg_K/s',  &
                        missing_value=missing_value               )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  & 
                        missing_value=missing_value               )

   id_rh_cmip = register_diag_field ( mod_name, &
     'rh_cmip', axes(1:3), Time, &
     'relative humidity',                            'percent',  &
                      missing_value=missing_value               )

   id_qs = register_diag_field ( mod_name, &
     'qs', axes(1:3), Time, &
         'saturation specific humidity',                 'kg/kg',    & 
                        missing_value=missing_value               )
   
if (do_donner_deep) then

   id_tdt_deep_donner= register_diag_field ( mod_name, &
           'tdt_deep_donner', axes(1:3), Time, &
           ' heating rate - deep portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_deep_donner = register_diag_field ( mod_name, &
           'qdt_deep_donner', axes(1:3), Time, &
           ' moistening rate - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_deep_donner = register_diag_field ( mod_name, &
     'qadt_deep_donner', axes(1:3), Time, &
     ' cloud amount tendency - deep portion', '1/s', &
                        missing_value=missing_value               )

   id_qldt_deep_donner = register_diag_field ( mod_name, &
     'qldt_deep_donner', axes(1:3), Time, &
     ' cloud liquid tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_deep_donner = register_diag_field ( mod_name, &
     'qidt_deep_donner', axes(1:3), Time, &
     ' ice water tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )
   if (do_liq_num) &
    id_qndt_deep_donner = register_diag_field ( mod_name, &
            'qndt_deep_donner', axes(1:3), Time, &
            'deep convection cloud drop tendency', '#/kg/s', &
                       missing_value=missing_value               )

   if (do_ice_num) &
     id_qnidt_deep_donner = register_diag_field ( mod_name, &
      'qnidt_deep_donner', axes(1:3), Time, &
     ' ice number tendency - deep portion', '#/kg/s', &
                         missing_value=missing_value               )

   id_tdt_mca_donner = register_diag_field ( mod_name, &
     'tdt_mca_donner', axes(1:3), Time, &
     ' heating rate - mca  portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_mca_donner = register_diag_field ( mod_name, &
           'qdt_mca_donner', axes(1:3), Time, &
           ' moistening rate - mca  portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_prec_deep_donner = register_diag_field ( mod_name, &
           'prc_deep_donner', axes(1:2), Time, &
           ' total precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_prec1_deep_donner = register_diag_field ( mod_name, &
           'prc1_deep_donner', axes(1:2), Time, &
           ' change in precip for conservation in donner', 'kg/m2/s ', &
              missing_value=missing_value, mask_variant = .true., &
             interp_method = "conserve_order1"  )

   id_prec_mca_donner = register_diag_field ( mod_name, &
           'prc_mca_donner', axes(1:2), Time, &
           ' total precip rate - mca  portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_deep_donner = register_diag_field ( mod_name, &
           'snow_deep_donner', axes(1:2), Time, &
           ' frozen precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_mca_donner = register_diag_field ( mod_name, &
           'snow_mca_donner', axes(1:2), Time, &
           ' frozen precip rate -  mca portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_mc_donner = register_diag_field ( mod_name, &
           'mc_donner', axes(1:3), Time, &
           'Net Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_mc_donner_half = register_diag_field ( mod_name, &
           'mc_donner_half', axes(half), Time, &
           'Net Mass Flux from donner at half levs',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_mc_conv_up = register_diag_field ( mod_name, &
           'mc_conv_up', axes(1:3), Time, &
          'Upward Mass Flux from convection',   'kg/m2/s', &
                       missing_value=missing_value               )

   id_m_cdet_donner = register_diag_field ( mod_name, &
           'm_cdet_donner', axes(1:3), Time, &
           'Detrained Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_m_cellup = register_diag_field ( mod_name, &
           'm_cellup', axes(half), Time, &
           'Upward Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

! ---> h1g, cell and meso-scale cloud fraction from donner deep, 2011-08-08
   id_cell_cld_frac = register_diag_field ( mod_name, &
           'cell_cld_frac', axes(1:3), Time, & 
           'cell cloud fraction from donner',   '', &
                        missing_value=missing_value               )

   id_meso_cld_frac = register_diag_field ( mod_name, &
           'meso_cld_frac', axes(1:3), Time, & 
           'meso-scale cloud fraction from donner',   '', &
                        missing_value=missing_value               )

   id_donner_humidity_area = register_diag_field ( mod_name, &
           'donner_humidity_area', axes(1:3), Time, &
           'donner humidity area',  '', &
                        missing_value=missing_value               )
! <--- h1g, cell and meso-scale cloud fraction from donner deep, 2011-08-08


endif


if (do_uw_conv) then

   id_tdt_uw = register_diag_field ( mod_name, &
           'tdt_uw', axes(1:3), Time, &
           'UW convection heating rate', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_uw = register_diag_field ( mod_name, &
           'qdt_uw', axes(1:3), Time, &
           'UW convection moistening rate', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_uw = register_diag_field ( mod_name, &
           'qadt_uw', axes(1:3), Time, &
           'UW convection cloud amount tendency', '1/s', &
                        missing_value=missing_value               )

   id_qldt_uw = register_diag_field ( mod_name, &
           'qldt_uw', axes(1:3), Time, &
           'UW convection cloud liquid tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_uw = register_diag_field ( mod_name, &
           'qidt_uw', axes(1:3), Time, &
           'UW convection ice water tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   if (do_liq_num) &
    id_qndt_uw = register_diag_field ( mod_name, &
           'qndt_uw', axes(1:3), Time, &
           'UW convection cloud drop tendency', '#/kg/s', &
                        missing_value=missing_value               )

    if (do_ice_num) &
     id_qnidt_uw = register_diag_field ( mod_name, &
           'qnidt_uw', axes(1:3), Time, &
           'UW convection ice number tendency', '#/kg/s', &
                        missing_value=missing_value               )

endif

   id_qvout = register_diag_field ( mod_name, &
           'qvout', axes(1:3), Time, 'qv after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qaout = register_diag_field ( mod_name, &
           'qaout', axes(1:3), Time, 'qa after strat_cloud', 'none', &
                        missing_value=missing_value               )

   id_qlout = register_diag_field ( mod_name, &
           'qlout', axes(1:3), Time, 'ql after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qiout = register_diag_field ( mod_name, &
           'qiout', axes(1:3), Time, 'qi after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   if (do_liq_num) then
   id_qnout = register_diag_field ( mod_name, &
           'qnout', axes(1:3), Time, 'qn after strat_cloud', '#/kg', &
                        missing_value=missing_value               )
   endif

   if (do_ice_num) then
   id_qniout = register_diag_field ( mod_name, &
           'qniout', axes(1:3), Time, 'qni after strat_cloud', '#/kg', &
                        missing_value=missing_value               )
   endif

!---------------------------------------------------------------------
!    register diagnostics for lightning NOx
!---------------------------------------------------------------------

   if (get_tracer_index(MODEL_ATMOS,'no') > 0) &
     id_prod_no = register_diag_field ( 'tracers', &
             'hook_no', axes(1:3), Time, &
             'hook_no',   'molec/cm3/s')

!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!    register the diagnostics associated with convective tracer 
!    transport.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv    (num_tracers))
      allocate (id_tracerdt_conv_col(num_tracers))
      allocate (id_wet_deposition(num_tracers))
      allocate (id_wetdep       (num_tracers))
      allocate (id_conv_tracer           (num_tracers))
      allocate (id_conv_tracer_col(num_tracers))

      id_tracerdt_conv = -1
      id_tracerdt_conv_col = -1
      id_wet_deposition = -1
      id_wetdep = -1
      id_conv_tracer = -1
      id_conv_tracer_col = -1
      
 
      id_wetdep_om = &
                         register_diag_field ( mod_name, &
                         'om_wet_dep',  &
                         axes(1:2), Time,  &
                         'total om wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_SOA = &
                         register_diag_field ( mod_name, &
                         'SOA_wet_dep',  &
                         axes(1:2), Time,  &
                         'total SOA wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_bc = &
                         register_diag_field ( mod_name, &
                         'bc_wet_dep',  &
                         axes(1:2), Time,  &
                         'total bc wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_so4 = &
                         register_diag_field ( mod_name, &
                         'so4_wet_dep',  &
                         axes(1:2), Time,  &
                         'total so4 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_so2 = &
                         register_diag_field ( mod_name, &
                         'so2_wet_dep',  &
                         axes(1:2), Time,  &
                         'total so2 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_DMS = &
                         register_diag_field ( mod_name, &
                         'DMS_wet_dep',  &
                         axes(1:2), Time,  &
                         'total DMS wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_NH4NO3 =  &
                         register_diag_field ( mod_name, &
                         'totNH4_wet_dep',  &
                         axes(1:2), Time,  &
                         'total NH4 + NH3 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_salt   =  &
                         register_diag_field ( mod_name, &
                         'ssalt_wet_dep',  &
                         axes(1:2), Time,  &
                         'total seasalt wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_dust   =  &
                         register_diag_field ( mod_name, &
                         'dust_wet_dep',  &
                         axes(1:2), Time,  &
                         'total dust wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_f_snow_berg   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg',  &
                         axes(1:3), Time,  &
                         'fraction of snow/ice produced having IFN', &
                         'fraction',  &
                         missing_value=missing_value)

      id_f_snow_berg_cond   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg_cond',  &
                         axes(1:3), Time,  &
                         'conditional fraction of snow/ice produced &
                         &having IFN', 'fraction',  &
                         mask_variant = .true., &
                         missing_value=missing_value)

      id_f_snow_berg_wtd   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg_wtd',  &
                         axes(1:3), Time,  &
                         'product of snow/ice produced having IFN and &
                         &ls precip falling out of gridbox', &
                         'kg(h2o)/m2/s', mask_variant = .true.,   &
                         missing_value=missing_value)

      do n = 1,num_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        if (tracers_in_donner(n) .or. &
            tracers_in_ras(n)      .or.  &
            tracers_in_mca(n)      .or.  &
            tracers_in_uw(n)) then

          diaglname = trim(tracer_name)//  &
                        ' wet deposition from all precip'
          id_wetdep(n) = &
                       register_diag_field ( mod_name, &
                          TRIM(tracer_name)//'_wet_depo',  &
                          axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'/s',  &
                          missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'*(kg/m2)/s',   &
                         missing_value=missing_value)
         endif
 
         diaglname = trim(tracer_name)
         id_conv_tracer(n) =    &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
         diaglname =  ' column integrated' // trim(tracer_name)
         id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
         id_wet_deposition(n) = register_diag_field( mod_name, &
           trim(tracer_name)//'_wetdep', axes(1:3), Time, &
           trim(tracer_name)//' tendency from wet deposition',TRIM(tracer_units)//'/sec', &
           missing_value=missing_value )
      end do

!------------------------------------------------------------------
!    register the variables associated with the mca component of 
!    donner_deep transport.
!------------------------------------------------------------------
     if (do_donner_deep) then
       allocate (id_tracerdt_mcadon  (num_donner_tracers))
       allocate (id_tracerdt_mcadon_col(num_donner_tracers))
 
       nn = 1
       do n = 1,num_tracers
         call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                units = tracer_units)
         if (tracers_in_donner(n) ) then
           diaglname = trim(tracer_name)//  &
                       ' tendency from donner-mca'
           id_tracerdt_mcadon(nn) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_donmca',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                        missing_value=missing_value)

           diaglname = trim(tracer_name)//  &
                       ' total path tendency from donner-mca'
           id_tracerdt_mcadon_col(nn) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_donmca_col', &
                        axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'*(kg/m2)/s',   &
                        missing_value=missing_value)
           nn = nn + 1
         endif
       end do
 

     endif

end subroutine diag_field_init


!#######################################################################
function doing_strat()
logical :: doing_strat

  if (.not. module_is_initialized) call error_mesg ('doing_strat',  &
                     'moist_processes_init has not been called.', FATAL)

  doing_strat = do_strat

end function doing_strat


!#######################################################################  

subroutine set_cosp_precip_sources (cosp_precip_sources)

character(len=16),        intent(in) :: cosp_precip_sources

     if (trim(cosp_precip_sources)  == 'stratdeepuw') then
       strat_precip_in_cosp = 1.
       donner_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'stratdeep') then
       strat_precip_in_cosp = 1.
       donner_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'stratuw') then
       strat_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'deepuw') then
       donner_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'strat') then
       strat_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'deep') then
       donner_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'uw') then
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'noprecip') then
!     COSP run without any precip input     
     else
       call error_mesg ('moist_processes_mod:set_cosp_precip_sources', &
        'cosp_precip_sources does not match any currently allowed string',&
                                                                 FATAL)
     endif

end subroutine set_cosp_precip_sources


!#######################################################################




end module moist_processes_mod

  

