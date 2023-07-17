module ml_interface

#ifdef INTERNAL_FILE_NML
    use mpp_mod, only: input_nml_file
#else
    use fms_mod, only: open_namelist_file, close_file
#endif

use fms_mod, only: write_version_number, file_exist, close_file, stdlog, error_mesg, NOTE, FATAL, read_data, field_size, uppercase, mpp_pe, check_nml_error

use interpolator_mod, only: interpolate_type,interpolator_init&
     &,CONSTANT,interpolator

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: ml_interface.F90,v 1.0'

character(len=128) :: tagname= &
'$Name:  $'
character(len=10), parameter :: mod_name='ml_interface'

!=================================================================================================================================

public :: ml_generated_file_init, read_ml_generated_file

character(len=256) :: conv_input_file  = 'ml_input'
character(len=256) :: tstd_field_name = 'tstd' 
character(len=256) :: qstd_field_name = 'qstd' 

namelist / ml_interface_nml / conv_input_file, tstd_field_name, qstd_field_name 


logical :: module_is_initialized =.false.
type(interpolate_type),save :: conv_input_file_interp


contains


subroutine ml_interface_init(is, ie, js, je, rad_lonb_2d, rad_latb_2d)
    
    real, intent(in), dimension(:,:) :: rad_lonb_2d, rad_latb_2d
    integer, intent(in) :: is, ie, js, je

    if(module_is_initialized) return

    call write_version_number(version, tagname)

    unit = open_namelist_file ()
    ierr=1
    do while (ierr /= 0)
    read  (unit, nml=ml_interface_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'ml_interface_nml')
    enddo
    10 call close_file (unit)

    if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=mixed_layer_nml)


    call interpolator_init( conv_input_file_interp, trim(conv_input_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )

    return
end subroutine ml_interface_init


subroutine read_ml_generated_file(Time, p_half, tstd, qstd)

    type(time_type), intent(in)         :: Time
    real, dimension(:,:,:), intent(in)  :: p_half
    real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)), intent(out)                   :: tstd, qstd

    if(.not.module_is_initialized) then
        call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
      endif

    call interpolator( conv_input_file_interp, p_half, tstd, tstd_field_name)
    call interpolator( conv_input_file_interp, p_half, qstd, qstd_field_name)


end subroutine mixed_layer
