module shallow_physics_mod

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! This program is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------

use               fms_mod, only: open_namelist_file,   &
                                 open_restart_file,    &
                                 file_exist,           &
                                 check_nml_error,      &
                                 error_mesg,           &
                                 FATAL, WARNING,       &
                                 write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 fms_init, fms_end,    &
                                 read_data,            &
                                 write_data,           &
                                 set_domain,           &
                                 close_file,           &
                                 stdlog

use         transforms_mod, only: get_sin_lat, get_cos_lat,  &
                                  get_deg_lon, get_deg_lat,  &
                                  get_wts_lat, &
                                  get_grid_domain, get_spec_domain, &
                                  grid_domain

use       time_manager_mod, only: time_type, get_time

use diag_manager_mod, only: diag_axis_init, register_static_field, register_diag_field, send_data

!========================================================================
implicit none
private
!========================================================================

public :: shallow_physics_init,    &
          shallow_physics,         &
          shallow_physics_end,     &
          phys_type


! version information 
!========================================================================
character(len=128) :: version = '$Id: shallow_physics.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
!========================================================================

type phys_type
   real, pointer, dimension(:,:)   :: empty=>NULL()
end type

logical :: module_is_initialized = .false.

integer :: is, ie, js, je

integer :: pe
logical :: root

real, allocatable, dimension(:) :: rad_lat, deg_lat, deg_lon, &
         sin_lat, cos_lat, wts_lat

real, allocatable, dimension(:,:) :: h_eq

real    :: kappa_m, kappa_t

real, allocatable, dimension(:) :: storm_time
integer :: storm_count = 0
real, allocatable, dimension(:) :: storm_lat, storm_lon

integer :: id_dt_hg_physical_forcing, id_dt_hg_rad_forcing

integer, allocatable, dimension(:)   :: seed        ! random number seed
integer :: nseed


! namelist 
!========================================================================

real    :: fric_damp_time  = -20.0
real    :: therm_damp_time = -10.0
real    :: del_h           = 1.e04
real    :: h_0             = 3.e04
real    :: h_amp           = 2.e04
real    :: h_lon           =  90.0
real    :: h_lat           =  25.0
real    :: h_width         =  15.0
real    :: h_itcz          = 1.e05
real    :: itcz_width      =  4.0

namelist /shallow_physics_nml/ fric_damp_time, therm_damp_time, del_h, h_0, &
                               h_amp, h_lon, h_lat, h_width, &
                               itcz_width, h_itcz
!========================================================================

contains

!========================================================================

subroutine shallow_physics_init(Phys, Time, id_lon, id_lat) 

type(phys_type), intent(inout) :: Phys
type(time_type), intent(in)    :: Time
integer, intent(in) :: id_lon, id_lat

integer :: i, j, unit, ierr, io

real :: xx, yy, dd

character(len=64) :: file

call write_version_number(version, tagname)

pe = mpp_pe()
root = (pe == mpp_root_pe())

! read the namelist

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr=1
  do while (ierr /= 0)
    read  (unit, nml=shallow_physics_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'shallow_physics_nml')
  enddo
  10 call close_file (unit)
endif

if(fric_damp_time  < 0.0)  fric_damp_time = -  fric_damp_time*86400
if(therm_damp_time < 0.0) therm_damp_time = - therm_damp_time*86400

kappa_m = 0.0
kappa_t = 0.0
if( fric_damp_time .ne. 0.0) kappa_m = 1./fric_damp_time
if(therm_damp_time .ne. 0.0) kappa_t = 1./therm_damp_time

call get_grid_domain(is,ie,js,je)

allocate ( rad_lat      (js:je) )
allocate ( deg_lat      (js:je) )
allocate ( sin_lat      (js:je) )
allocate ( cos_lat      (js:je) )
allocate ( wts_lat      (js:je) )
allocate ( deg_lon      (is:ie) )
allocate ( h_eq   (is:ie,js:je) )

call get_wts_lat(wts_lat)
call get_deg_lat(deg_lat)
call get_deg_lon(deg_lon)
rad_lat = deg_lat*atan(1.)/45. 
sin_lat = sin(rad_lat)
cos_lat = cos(rad_lat)


do j = js, je
  do i = is, ie
     xx = (deg_lon(i) - h_lon)/(h_width*2.0)
     yy = (deg_lat(j) - h_lat)/h_width
     dd =  xx*xx + yy*yy
     h_eq(i,j) = h_0 + h_amp*max(1.e-10, exp(-dd))
  end do
end do

do j = js, je
  yy = deg_lat(j)/itcz_width
  dd = yy*yy
  h_eq(:,j) = h_eq(:,j) + h_itcz*exp(-dd)
end do 

id_dt_hg_physical_forcing       = register_diag_field  ('physics_mod', 'dt_hg_physical_forcing',     (/id_lon,id_lat/), Time, 'dt_hg_physical_forcing', 'm/sec')
id_dt_hg_rad_forcing       = register_diag_field  ('physics_mod', 'dt_hg_rad_forcing',     (/id_lon,id_lat/), Time, 'dt_hg_rad_forcing', 'm/sec')

call random_seed(size=nseed)
allocate(seed(nseed))

allocate ( storm_time  (0:30)) ; storm_time = 0.
allocate ( storm_lat   (0:30)) ; storm_lat  = 0.
allocate ( storm_lon   (0:30)) ; storm_lon  = 0.

if(file_exist('INPUT/shallow_physics.res')) then
 file = 'INPUT/shallow_physics.res'
 call read_data(trim(file), 'storm_time', storm_time)
 call read_data(trim(file), 'storm_lat', storm_lat)
 call read_data(trim(file), 'storm_lon', storm_lon) 
 call read_data(trim(file), 'ran_nmbr_seed_ph_forcing', seed, no_domain=.true.)
 call random_seed(put=seed) 
endif

module_is_initialized = .true.

return
end subroutine shallow_physics_init

!=======================================================================

subroutine shallow_physics(Time, dt_ug, dt_vg, dt_hg, ug, vg, hg,   &
                             delta_t, previous, current, Phys)

real, intent(inout),  dimension(is:ie, js:je)    :: dt_ug, dt_vg, dt_hg
real, intent(in)   ,  dimension(is:ie, js:je, 2) :: ug, vg, hg

real   , intent(in)  :: delta_t
integer, intent(in)  :: previous, current

type(time_type), intent(in)    :: Time
type(phys_type), intent(inout) :: Phys

integer :: seconds,days
integer :: i, j, unit, ierr, io, ii, jj
    real :: xx, yy, dd, deg_lon0, deg_lat0, rad_lon0, rad_lat0, mm, tt
    real :: storm_length
    real :: storm_interval
    real :: storm_strength
    real :: h_width =  2.0
    real :: model_time, gs, te, ke, pe
    ! real, dimension(0:30), save :: storm_time = 0
    ! integer, save :: storm_count = 0
    integer :: storm_count_i = 0
    ! real, dimension(0:30), save :: storm_lat, storm_lon

logical :: used

real, dimension(is:ie, js:je) :: dt_hg_physical_forcing, dt_hg_rad_forcing

dt_hg_physical_forcing = 0.

call get_time(Time,seconds,days)

storm_strength=1.0
model_time = days*86400+seconds

storm_interval = 100000.0!0.5*(10**therm_damp_time) / 100.0
storm_length = 100000.0!0.5*(10**therm_damp_time) / 100.0


!call random_seed()
if (model_time == 0.0) then
        storm_count = 0
        storm_strength = 0.0
        storm_time(0) = storm_interval * 1.5
else if (mod(model_time, storm_interval) == 0) then
            if (storm_count == 30) then
                    storm_count = 0
            else
                storm_count = storm_count + 1
            endif    

            if (storm_count == 0) then
                storm_time(storm_count) = storm_time(30) + storm_interval
            else
                storm_time(storm_count) = storm_time(storm_count-1) + storm_interval
            endif 

            call random_number(storm_lon(storm_count))
            call random_number(storm_lat(storm_count))
            storm_lon(storm_count) = storm_lon(storm_count)* 360.
            storm_lat(storm_count) = - (90. - 45.*acos(2*storm_lat(storm_count)-1)/atan(1.))
end if

do storm_count_i = 0,30
tt = ((model_time - storm_time(storm_count_i))**2)/storm_length**2
storm_strength =   1.0 * (h_eq(is,js)) / storm_length
call get_wts_lat(wts_lat)
call get_deg_lat(deg_lat)
call get_deg_lon(deg_lon)
rad_lat = deg_lat*atan(1.)/45.
sin_lat = sin(rad_lat)
cos_lat = cos(rad_lat)
      do mm = 0, 1
      do j = js, je
         do i = is, ie
            xx = (deg_lon(i) - (storm_lon(storm_count_i)+mm*360.))/(h_width/cos_lat(j))
            yy = (deg_lat(j) - storm_lat(storm_count_i))/h_width
            dd =  xx*xx + yy*yy
            if (dd < 4 * h_width) then
               dt_hg_physical_forcing(i,j) = storm_strength * exp(-dd) * exp(-tt)
               dt_hg(i,j) = dt_hg(i,j) + dt_hg_physical_forcing(i,j)
            end if
         end do
      end do
      end do
end do


dt_ug = dt_ug - kappa_m*ug(:,:,previous)
dt_vg = dt_vg - kappa_m*vg(:,:,previous)

dt_hg_rad_forcing = - kappa_t*(hg(:,:,previous) - h_eq)
dt_hg = dt_hg + dt_hg_rad_forcing

used = send_data(id_dt_hg_physical_forcing, dt_hg_physical_forcing, Time)
used = send_data(id_dt_hg_rad_forcing, dt_hg_rad_forcing, Time)

return
end subroutine shallow_physics

!======================================================================

subroutine shallow_physics_end(Phys)

type(phys_type), intent(in) :: Phys

integer :: unit

character(len=64) :: file

if(.not.module_is_initialized) then
  call error_mesg('shallow_physics_end','physics has not been initialized ', FATAL)
endif

file='RESTART/shallow_physics.res'
call write_data(trim(file), 'storm_time', storm_time)
call write_data(trim(file), 'storm_lat', storm_lat)
call write_data(trim(file), 'storm_lon', storm_lon)
call random_seed(get=seed)
call write_data(trim(file), 'ran_nmbr_seed_ph_forcing', seed, no_domain=.true.)
module_is_initialized = .false.

return
end subroutine shallow_physics_end

!======================================================================

end module shallow_physics_mod
