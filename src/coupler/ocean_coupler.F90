!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of Isca.                                        *
!*                                                                     *
!* Isca is free software; you can redistribute it and/or modify it     *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!***********************************************************************

module ocean_coupler_mod

! Exchange of surface fields between Isca and an external ocean component
! (MOM6, an ML ocean emulator, or a test slab ocean) through OASIS3-MCT.
!
! Isca and the ocean run as separate executables launched together
! (mpirun -np N isca.x : -np M ocean.x). OASIS handles the exchange,
! time-averaging (LOCTRANS AVERAGE) and remapping between the grids, as
! configured in the 'namcouple' file in the run directory.
!
! Everything that needs the OASIS library is inside #ifdef ISCA_OASIS, so
! Isca builds without OASIS are unaffected. Asking for coupling in a build
! without ISCA_OASIS is a FATAL error.
!
! Sign conventions of fields sent to the ocean (all on the Isca grid):
!   ISCA_TAUX, ISCA_TAUY : wind stress ON THE OCEAN [Pa], positive eastward/northward
!   ISCA_SWNET           : net shortwave into the ocean [W/m2], positive down
!   ISCA_LWNET           : net longwave into the ocean [W/m2], positive down
!   ISCA_SENS            : sensible heat flux [W/m2], positive UPWARD (out of ocean)
!   ISCA_EVAP            : evaporation [kg/m2/s], positive UPWARD (out of ocean)
!   ISCA_PRECIP          : liquid precipitation [kg/m2/s], positive down
!   ISCA_PSURF           : surface pressure [Pa]
!   ISCA_RUNOFF          : river runoff [kg/m2/s], positive into the ocean. Water that
!                          leaves the surface anywhere the ocean is not (bucket overflow
!                          on land, P-E over uncoupled slab points) is moved to an ocean
!                          cell given by INPUT/runoff_routing.nc, so no water is lost.
! Fields received from the ocean:
!   OCN_SST              : sea-surface temperature [K]
!   OCN_FRAC             : fraction of the Isca cell covered by coupled ocean [0-1]
!   OCN_FRAZIL           : (optional) heat the ocean added to stop water cooling below
!                          freezing [W/m2 of ocean], averaged over the last coupling period.
!                          With Isca's thermodynamic sea ice this heat grows ice instead.

#ifdef ISCA_OASIS
use mod_oasis
#endif

use            fms_mod, only: error_mesg, FATAL, NOTE, WARNING, mpp_pe, mpp_root_pe, &
                              file_exist, write_version_number
use            mpp_mod, only: mpp_sum
use   time_manager_mod, only: time_type, get_time, operator(-)
use      constants_mod, only: RADIUS, PI

implicit none
private

public :: ocean_coupler_init_comp, ocean_coupler_init, ocean_coupler_put, &
          ocean_coupler_get, ocean_coupler_end, ocean_coupler_terminate,  &
          ocean_coupler_active, ocean_coupler_put_runoff

character(len=128), parameter :: version = '$Id: ocean_coupler.F90 2026-10-23 $'
character(len=128), parameter :: tagname = '$Name: coupled_ocean_oasis $'

integer, parameter :: n_out = 9, n_in = 3
character(len=11), parameter :: out_names(n_out) = (/ 'ISCA_TAUX  ', 'ISCA_TAUY  ', &
     'ISCA_SWNET ', 'ISCA_LWNET ', 'ISCA_SENS  ', 'ISCA_EVAP  ', 'ISCA_PRECIP', 'ISCA_PSURF ', &
     'ISCA_RUNOFF' /)
character(len=10), parameter :: in_names(n_in) = (/ 'OCN_SST   ', 'OCN_FRAC  ', 'OCN_FRAZIL' /)
character(len=4),  parameter :: grid_name = 'isca'
character(len=4),  parameter :: grid_name_runoff = 'iscr'   ! unmasked copy, for runoff
real, parameter :: frac_min_routed = 1.e-3   ! least ocean fraction a cell needs to receive routed runoff

logical :: oasis_active = .false.   ! OASIS was initialised (a namcouple exists)
logical :: module_is_initialized = .false.
integer :: comp_id = -1, part_id = -1
integer :: out_ids(n_out) = -1, in_ids(n_in) = -1
integer :: is, ie, js, je, nlon_g, nlat_g
type(time_type) :: Time_run_start
real, allocatable :: cell_area(:,:)          ! local cell areas [m2]
integer, allocatable :: runoff_dest(:,:)     ! global index (i + (j-1)*nlon) a cell's runoff goes to; 0 = none
real, allocatable :: frac_last(:,:)          ! ocean fraction of each cell, as last received

contains

!#######################################################################

subroutine ocean_coupler_init_comp(localcomm, coupled)
! Must be called before fms_init. If a 'namcouple' file is present, start
! OASIS and return the communicator Isca must use. Otherwise return
! coupled = .false. and leave MPI initialisation to FMS as usual.
integer, intent(out) :: localcomm
logical, intent(out) :: coupled
#ifdef ISCA_OASIS
integer :: ierr
logical :: have_namcouple

inquire(file='namcouple', exist=have_namcouple)
coupled = have_namcouple
localcomm = -1
if (.not. have_namcouple) return

call oasis_init_comp(comp_id, 'isca', ierr)
if (ierr /= 0) call oasis_abort(comp_id, 'ocean_coupler_init_comp', 'oasis_init_comp failed')
call oasis_get_localcomm(localcomm, ierr)
if (ierr /= 0) call oasis_abort(comp_id, 'ocean_coupler_init_comp', 'oasis_get_localcomm failed')
oasis_active = .true.
#else
coupled = .false.
localcomm = -1
#endif
end subroutine ocean_coupler_init_comp

!#######################################################################

logical function ocean_coupler_active()
ocean_coupler_active = oasis_active
end function ocean_coupler_active

!#######################################################################

subroutine ocean_coupler_init(Time, is_in, ie_in, js_in, je_in, num_lon, num_lat, &
                              rad_lonb_2d, rad_latb_2d, land)
! Define the Isca grid and partition to OASIS, declare the coupling fields
! and finish the OASIS definition phase (collective with the ocean).
type(time_type), intent(in) :: Time
integer, intent(in) :: is_in, ie_in, js_in, je_in, num_lon, num_lat
real,    intent(in) :: rad_lonb_2d(:,:), rad_latb_2d(:,:)   ! local cell corners [radians]
logical, intent(in) :: land(:,:)
#ifdef ISCA_OASIS
real, dimension(num_lon+1) :: lonb_g, cnt_lon
real, dimension(num_lat+1) :: latb_g, cnt_lat
real(kind=8), allocatable :: lon_c(:,:), lat_c(:,:), clon(:,:,:), clat(:,:,:), area(:,:)
integer, allocatable :: mask(:,:)
integer :: ig_paral(5), var_nodims(2), var_shape(4), ierr, n, i, j, iflag
real, dimension(num_lon, num_lat) :: land_g

if (module_is_initialized) return
if (.not. oasis_active) call error_mesg('ocean_coupler_init', &
    'ocean coupling requested but no namcouple file was found in the run directory', FATAL)
call write_version_number(version, tagname)

is = is_in; ie = ie_in; js = js_in; je = je_in
nlon_g = num_lon; nlat_g = num_lat
Time_run_start = Time

! Assemble the global cell boundaries from the local ones. Boundaries shared
! between processors are counted more than once, so divide by the count.
lonb_g = 0.; cnt_lon = 0.; latb_g = 0.; cnt_lat = 0.
lonb_g(is:ie+1) = rad_lonb_2d(:,1);  cnt_lon(is:ie+1) = 1.
latb_g(js:je+1) = rad_latb_2d(1,:);  cnt_lat(js:je+1) = 1.
call mpp_sum(lonb_g, num_lon+1); call mpp_sum(cnt_lon, num_lon+1)
call mpp_sum(latb_g, num_lat+1); call mpp_sum(cnt_lat, num_lat+1)
lonb_g = lonb_g / cnt_lon
latb_g = latb_g / cnt_lat

allocate(cell_area(is:ie, js:je), frac_last(is:ie, js:je))
frac_last = 1.
do j = js, je
  do i = is, ie
    cell_area(i,j) = RADIUS**2 * (lonb_g(i+1) - lonb_g(i)) * (sin(latb_g(j+1)) - sin(latb_g(j)))
  enddo
enddo
call read_runoff_routing()

land_g = 0.
where (land) land_g(is:ie, js:je) = 1.
call mpp_sum(land_g, num_lon*num_lat)

! Grid description, written once by the root PE (grids.nc, masks.nc, areas.nc).
! The Isca mask marks land as inactive; for an aquaplanet everything is active.
if (mpp_pe() == mpp_root_pe()) then
  allocate(lon_c(num_lon,num_lat), lat_c(num_lon,num_lat), area(num_lon,num_lat), &
           clon(num_lon,num_lat,4), clat(num_lon,num_lat,4), mask(num_lon,num_lat))
  do j = 1, num_lat
    do i = 1, num_lon
      lon_c(i,j) = 0.5*(lonb_g(i) + lonb_g(i+1)) * 180./PI
      lat_c(i,j) = 0.5*(latb_g(j) + latb_g(j+1)) * 180./PI
      ! corners counter-clockwise from south-west
      clon(i,j,:) = (/ lonb_g(i), lonb_g(i+1), lonb_g(i+1), lonb_g(i) /) * 180./PI
      clat(i,j,:) = (/ latb_g(j), latb_g(j),   latb_g(j+1), latb_g(j+1) /) * 180./PI
      area(i,j) = RADIUS**2 * (lonb_g(i+1) - lonb_g(i)) * (sin(latb_g(j+1)) - sin(latb_g(j)))
    enddo
  enddo
  mask = nint(land_g)
  call oasis_start_grids_writing(iflag)
  call oasis_write_grid(grid_name, num_lon, num_lat, lon_c, lat_c)
  call oasis_write_corner(grid_name, num_lon, num_lat, 4, clon, clat)
  call oasis_write_mask(grid_name, num_lon, num_lat, mask)
  call oasis_write_area(grid_name, num_lon, num_lat, area)
  ! The same grid with no mask, for runoff: remapped with DESTAREA normalisation it
  ! delivers to the ocean exactly the ocean-covered fraction of each cell's value
  mask = 0
  call oasis_write_grid(grid_name_runoff, num_lon, num_lat, lon_c, lat_c)
  call oasis_write_corner(grid_name_runoff, num_lon, num_lat, 4, clon, clat)
  call oasis_write_mask(grid_name_runoff, num_lon, num_lat, mask)
  call oasis_write_area(grid_name_runoff, num_lon, num_lat, area)
  call oasis_terminate_grids_writing()
  deallocate(lon_c, lat_c, area, clon, clat, mask)
endif

! Box partition: this PE owns i=is:ie, j=js:je of the global grid.
ig_paral(1) = 2
ig_paral(2) = (js-1)*num_lon + (is-1)
ig_paral(3) = ie - is + 1
ig_paral(4) = je - js + 1
ig_paral(5) = num_lon
call oasis_def_partition(part_id, ig_paral, ierr, num_lon*num_lat)

var_nodims = (/ 2, 1 /)
var_shape  = (/ 1, ie-is+1, 1, je-js+1 /)
do n = 1, n_out
  call oasis_def_var(out_ids(n), trim(out_names(n)), part_id, var_nodims, OASIS_Out, &
                     var_shape, OASIS_Real, ierr)
enddo
do n = 1, n_in
  call oasis_def_var(in_ids(n), trim(in_names(n)), part_id, var_nodims, OASIS_In, &
                     var_shape, OASIS_Real, ierr)
enddo
if (in_ids(1) < 0) call error_mesg('ocean_coupler_init', &
    'OCN_SST is not in the namcouple: the ocean SST cannot be received', FATAL)

call oasis_enddef(ierr)
if (ierr /= 0) call oasis_abort(comp_id, 'ocean_coupler_init', 'oasis_enddef failed')

module_is_initialized = .true.
#else
call error_mesg('ocean_coupler_init', 'ocean coupling requested, but Isca was compiled '// &
     'without -DISCA_OASIS (no OASIS3-MCT support)', FATAL)
#endif
end subroutine ocean_coupler_init

!#######################################################################

subroutine ocean_coupler_put(Time, taux, tauy, swnet, lwnet, sens, evap, precip, psurf)
! Called every atmospheric timestep. OASIS averages the puts over the
! coupling period and sends the average when date+LAG is a coupling time.
type(time_type), intent(in) :: Time
real, intent(in), dimension(:,:) :: taux, tauy, swnet, lwnet, sens, evap, precip, psurf
#ifdef ISCA_OASIS
integer :: date, info

date = seconds_since_start(Time)
call put_one(1, taux);  call put_one(2, tauy);  call put_one(3, swnet)
call put_one(4, lwnet); call put_one(5, sens);  call put_one(6, evap)
call put_one(7, precip); call put_one(8, psurf)

contains

subroutine put_one(n, field)
integer, intent(in) :: n
real, intent(in) :: field(:,:)
real(kind=8) :: buf(size(field,1), size(field,2))
if (out_ids(n) < 0) return   ! field not in the namcouple
buf = field
call oasis_put(out_ids(n), date, buf, info)
end subroutine put_one
#endif
end subroutine ocean_coupler_put

!#######################################################################

subroutine ocean_coupler_put_runoff(Time, runoff, routed, spread)
! Route runoff [kg/m2/s] from the cells where it is produced to their
! destination ocean cells and send it. Called every atmospheric timestep.
! Routing is a global operation (rivers cross PE boundaries), done with one
! mpp_sum of a global array: cheap at the resolutions Isca is used at.
! spread [kg/m2/s]: water (of either sign) that is not river runoff, e.g. P-E over
! parts of cells the coupled ocean does not cover. It is added uniformly over the
! whole coupled ocean: routing it would concentrate it, negative when evaporation
! wins, in single coastal cells (salinifying and draining a shallow river mouth).
type(time_type), intent(in) :: Time
real, intent(in), dimension(:,:) :: runoff
real, intent(out), dimension(:,:), optional :: routed   ! runoff arriving in each cell [kg/m2/s]
real, intent(in), dimension(:,:), optional :: spread
#ifdef ISCA_OASIS
real :: routed_g(nlon_g*nlat_g), sums(2)
real(kind=8) :: buf(is:ie, js:je)
integer :: i, j, k, date, info

if (present(routed)) routed = 0.
if (out_ids(9) < 0) return
routed_g = 0.
do j = js, je
  do i = is, ie
    ! dest 0 marks an ocean cell: its own runoff (e.g. P-E over prescribed sea ice,
    ! which is not coupled at that time) stays in the cell
    k = runoff_dest(i,j)
    if (k == 0) k = i + (j-1)*nlon_g
    routed_g(k) = routed_g(k) + runoff(i-is+1, j-js+1) * cell_area(i,j)   ! [kg/s]
  enddo
enddo
call mpp_sum(routed_g, nlon_g*nlat_g)
! total water to spread [kg/s] and the coupled ocean area [m2]
sums = 0.
if (present(spread)) sums(1) = sum(spread * cell_area)
sums(2) = sum(min(max(frac_last, 0.), 1.) * cell_area)
do j = js, je
  do i = is, ie
    k = i + (j-1)*nlon_g
    if (frac_last(i,j) < frac_min_routed) then
      ! a destination the ocean does not cover (the routing map follows Isca's land
      ! mask): its runoff would be lost, so spread it with the rest
      sums(1) = sums(1) + routed_g(k)
      buf(i,j) = 0.
    else
      ! Remapping (unmasked grid, DESTAREA) delivers only the ocean-covered fraction of
      ! a cell's value to the ocean, so spread the routed water over that fraction to
      ! deliver all of it.
      buf(i,j) = routed_g(k) / (cell_area(i,j) * frac_last(i,j))
    endif
  enddo
enddo
call mpp_sum(sums, 2)
if (sums(2) > 0.) buf = buf + sums(1) / sums(2)
if (present(routed)) routed = buf * frac_last   ! per unit cell area
date = seconds_since_start(Time)
call oasis_put(out_ids(9), date, buf, info)
#else
if (present(routed)) routed = 0.
#endif
end subroutine ocean_coupler_put_runoff

!#######################################################################

subroutine ocean_coupler_get(Time, sst, ocean_frac, received, frazil)
! Called every atmospheric timestep. sst, ocean_frac and frazil are only changed
! when a coupling field actually arrives (received = .true.).
type(time_type), intent(in) :: Time
real, intent(inout), dimension(:,:) :: sst, ocean_frac
logical, intent(out) :: received
real, intent(inout), dimension(:,:), optional :: frazil   ! [W/m2 of ocean]; unchanged if the ocean sends none
#ifdef ISCA_OASIS
real(kind=8) :: buf(size(sst,1), size(sst,2))
integer :: date, info

received = .false.
date = seconds_since_start(Time)

call oasis_get(in_ids(1), date, buf, info)
if (is_received(info)) then
  sst = buf
  received = .true.
endif
if (in_ids(2) >= 0) then
  call oasis_get(in_ids(2), date, buf, info)
  if (is_received(info)) then
    ocean_frac = buf
    frac_last = buf
  endif
endif
if (present(frazil) .and. in_ids(3) >= 0) then
  call oasis_get(in_ids(3), date, buf, info)
  if (is_received(info)) frazil = buf
endif
#else
received = .false.
#endif
end subroutine ocean_coupler_get

!#######################################################################

subroutine ocean_coupler_end()
#ifdef ISCA_OASIS
module_is_initialized = .false.
#endif
end subroutine ocean_coupler_end

!#######################################################################

subroutine ocean_coupler_terminate()
! Called after fms_end: oasis_terminate finalises MPI.
#ifdef ISCA_OASIS
integer :: ierr
if (oasis_active) call oasis_terminate(ierr)
#endif
end subroutine ocean_coupler_terminate

!#######################################################################

#ifdef ISCA_OASIS
subroutine read_runoff_routing()
! INPUT/runoff_routing.nc: integer variable 'dest' (lon, lat), the global
! index i + (j-1)*nlon of the ocean cell each cell's runoff is sent to
! (0: an ocean cell, which keeps any runoff it produces itself). Written by
! isca.coupling. Without it, runoff is not routed (and not sent).
use netcdf
integer :: ncid, vid, status
integer, allocatable :: dest_g(:,:)

allocate(runoff_dest(is:ie, js:je))
runoff_dest = 0
if (.not. file_exist('INPUT/runoff_routing.nc')) then
  call error_mesg('ocean_coupler', 'no INPUT/runoff_routing.nc: runoff will not be routed '// &
                  'to the ocean (water over land/uncoupled points is lost)', WARNING)
  return
endif
allocate(dest_g(nlon_g, nlat_g))
status = nf90_open('INPUT/runoff_routing.nc', NF90_NOWRITE, ncid)
if (status == NF90_NOERR) status = nf90_inq_varid(ncid, 'dest', vid)
if (status == NF90_NOERR) status = nf90_get_var(ncid, vid, dest_g)
if (status /= NF90_NOERR) call error_mesg('ocean_coupler', 'reading INPUT/runoff_routing.nc: '// &
                                          trim(nf90_strerror(status)), FATAL)
status = nf90_close(ncid)
runoff_dest = dest_g(is:ie, js:je)
deallocate(dest_g)
end subroutine read_runoff_routing

integer function seconds_since_start(Time)
type(time_type), intent(in) :: Time
integer :: days, seconds
call get_time(Time - Time_run_start, seconds, days)
seconds_since_start = days*86400 + seconds
end function seconds_since_start

logical function is_received(info)
integer, intent(in) :: info
is_received = (info == OASIS_Recvd .or. info == OASIS_FromRest .or. &
               info == OASIS_RecvOut .or. info == OASIS_FromRestOut .or. &
               info == OASIS_Input)
end function is_received
#endif

end module ocean_coupler_mod
