! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE llcs
  USE qsat_wat_data, ONLY: T_low, T_high, delta_T, zerodegc, es

! Description:
!  Module containing all subroutines for the Lambert-Lewis
!  Convection Scheme (LLCS)
!  Further description of the scheme can be found below

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

!USE planet_constants_mod, ONLY: cp,r,kappa,rv,g,repsilon,               &
!     one_minus_epsilon, p_zero
!USE timestep_mod, ONLY: timestep
!USE water_constants_mod, ONLY: lc,hcapv

IMPLICIT NONE


! The only public interface is the 'UM Wrapper' module llcscontrol.
!PUBLIC llcs_control, calc_temp, cp, nlevels, r, g, timestep, lc

! -- PARAMETERS -- !
! Parameters (of the scheme / model)
REAL, PARAMETER :: timescale = 3600.0 ! Convective "timescale" in seconds
                                     ! (Tompkins and Craig (1998))

REAL, PARAMETER :: rh_crit = 0.8 ! Critical relative humidity required
! to trigger moist convection.

REAL, PARAMETER :: repsilon = 0.61298
REAL, PARAMETER :: one_minus_epsilon = 1 - repsilon
REAL, PARAMETER :: cp = 1005.0
REAL, PARAMETER :: r = 287.05
REAL, PARAMETER :: rv = r / repsilon
REAL, PARAMETER :: kappa = r / cp
REAL, PARAMETER :: g = 9.80655
REAL, PARAMETER :: lc = 2.501e6
REAL, PARAMETER :: hcapv = 1850.0
REAL, PARAMETER :: timestep = 1200.0
REAL, PARAMETER :: p_zero = 100000.0

!CHARACTER(LEN=*), PARAMETER :: ModuleName = 'LLCS'

CONTAINS

! ----------------------------------------------------------------------
! -- UM WRAPPER -- !
! ----------------------------------------------------------------------

SUBROUTINE llcs_control(t_s_array, q_s_array, p_surf_array,             &
p_half_array, p_full_array, t_f_array, q_f_array, rain_array)

!USE yomhook, ONLY: lhook, dr_hook
!USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! INPUTS
REAL, INTENT(IN) :: t_s_array(:,:,:), q_s_array(:,:,:),                 &
                    p_half_array(:,:,:), p_full_array(:,:,:)
REAL, INTENT(IN) :: p_surf_array(:,:)

! USED TO CALL CONVECTION SCHEME
INTEGER :: a, b, dim1, dim2, nlevels

! OUTPUTS
REAL, INTENT(OUT) :: t_f_array(:,:,:), q_f_array(:,:,:)
REAL, INTENT(OUT) :: rain_array(:,:)

!CHARACTER(LEN=*), PARAMETER :: RoutineName = 'LLCS_CONTROL'

!INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
!INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!REAL(KIND=jprb)               :: zhook_handle

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nlevels = SIZE(t_s_array,3)
dim1 = SIZE(t_s_array,1)
dim2 = SIZE(t_s_array,2)

DO a = 1 , dim1
  DO b = 1 , dim2

    CALL convection_scheme(t_s_array(a,b,:), q_s_array(a,b,:),          &
    p_surf_array(a,b), p_half_array(a,b,:), p_full_array(a,b,:),        &
    nlevels, t_f_array(a,b,:), q_f_array(a,b,:), rain_array(a,b))
  END DO
END DO

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE llcs_control

! ----------------------------------------------------------------------
! -- THE CONVECTION SCHEME -- !
! ----------------------------------------------------------------------

SUBROUTINE convection_scheme(theta_start, q_start, p_surf, p_half,      &
p_full, nlevels, theta_final, q_final, rain_sec)

! Required inputs are potential temperature, specific humidity, pressure
! on full model levels, pressure on half model levels and surface
! pressure.
! Potential temperature and specific humidity have (nlevels) levels.
! Exner pressure has (nlevels+1) levels.
! The arrays are set up such that, for example, theta_start(1) is
! potential temperature at the lowest model level and theta_start(nlevels) is
! potential temperature at the top of atmosphere.

! The convection scheme will output final potential temperature and
! final specific humidity, along with rainfall / second (kgm^-2s^-1)
! within the column.

! VARIABLE DESCRIPTIONS
!
! INPUT VARIABLES
! theta_start = inital potential temperature at full model levels.
! q_start = initial specific humidity at full model levels.
! p_full = pressure at full model levels
! p_half = pressure at half model levels
! nlevels = number of full model levels.
!
! CALCULATION VARIABLES
! theta_adjust = theta after correction to reference profile
! theta_noq_adjust = dummy theta profile where only dry correction takes
! place, used for enthalpy conservation.
! theta_adj_rlx = theta_adjust corrected to account for timestep /
! timescale.
! theta_noq_adj_rlx = as above.
! q_adjust = specific humidity after correction to q_sat.
! q_adj_rlx = q_adjust corrected to account for timestep / timescale.
! q_conserve = conserves specific humidity within the column (other than
! moisture lost through rainfall)
! q_dummy = dummy q profile to ensure moisture is not mixed twice
! convbase = the lowest level from which any form of mixing occurs.
! cloudbase  = a level from which adjustment to a moist adiabat can occur
! final_level = the final level of a convective event. Set to -1 whilst
! the event is in progress
! moist_has_happened = flag indicating whether moist convection has
! taken place within this column. -1 for no, 1 for yes.
! convection_flag = flag indicating whether convection has happened. It is
! initialised to -1 (no convection), and set to 1 if convection happens.
! e_sat = saturation vapour pressure.
! q_sat = saturation specific humidity.
!
! OUTPUT VARIABLES
! theta_final = final potential temperature at full model levels.
! q_final = final specific humidity at full model levels
! rainsec = convective precipitation within the column (kgm^-2s^-1)
!
! END OF VARIABLE DESCRIPTIONS

IMPLICIT NONE

! INPUTS
INTEGER, INTENT(IN) :: nlevels
REAL, INTENT(IN) :: theta_start(:), q_start(:), p_full(:), p_half(:)
REAL, INTENT(IN) :: p_surf

! FOR CALCULATION
REAL :: theta_adjust(nlevels), theta_noq_adjust(nlevels),               &
        theta_adj_rlx(nlevels), theta_noq_adj_rlx(nlevels),             &
        q_adjust(nlevels), q_adj_rlx(nlevels), q_conserve(nlevels),     &
        q_dummy(nlevels)
INTEGER :: convbase, cloudbase, final_level, j, moist_has_happened,     &
     convection_flag, moist_trigger
REAL :: rh_use

! OUTPUTS
REAL, INTENT(OUT) :: theta_final(:), q_final(:)
REAL, INTENT(OUT) :: rain_sec

! The scheme calls subroutines below here

! 1, Initialise variables.
CALL initialise(theta_start, theta_adjust, theta_noq_adjust,            &
theta_adj_rlx, theta_noq_adj_rlx, theta_final, q_start, q_adjust,       &
q_adj_rlx, q_conserve, q_dummy, q_final, rain_sec,                      &
final_level, moist_has_happened, convection_flag)

! 2 & 3), Diagnose and then perform convection.
DO j = 1, (nlevels-1)
  IF (final_level  ==  j) THEN
    CALL diagnose(nlevels, j, theta_adjust, q_adjust, p_full,           &
         p_surf, convbase, cloudbase, final_level, convection_flag,     &
         moist_trigger, rh_use)
    IF (convbase  /=  -1) THEN
      CALL convection(nlevels, theta_adjust, theta_noq_adjust, q_adjust,&
      q_dummy, p_full, p_surf, convbase, cloudbase, final_level,        &
      moist_has_happened, moist_trigger, rh_use)
    END IF
    ! 4 & 5), Dry and moist convection will be administered by the
    ! routine convection (above).
  END IF
END DO

IF (convection_flag  ==  1) THEN ! convection has happened so we need to
                                 ! perform the following steps.

  ! 6), Timescale
  ! If moist convection happens then we need to perform this step twice,
  ! once for the actual theta profile, and once for the 'dry convection
  ! only' profile used for enthalpy conservation
  IF (moist_has_happened  ==  1) THEN
    CALL relax(theta_start, theta_adjust, theta_adj_rlx, q_start,       &
    q_adjust, q_adj_rlx)
    CALL relax(theta_start, theta_noq_adjust, theta_noq_adj_rlx,        &
    q_start, q_adjust, q_dummy)
  ELSE
    CALL relax(theta_start, theta_adjust, theta_adj_rlx, q_start,       &
    q_adjust, q_adj_rlx)
  END IF

  ! 7), Conserve Moisture
  CALL moisture_conserve(nlevels, q_start, q_adj_rlx, p_full,           &
  p_half, p_surf, q_conserve)

  ! 8), Rain / Latent Heating
  IF (moist_has_happened  ==  1) THEN ! we need to consider rainfall from
                                      ! latent heating.

    CALL rain(nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx,   &
         q_conserve, p_full, p_half, p_surf, q_final, rain_sec)

  END IF

  ! 9), Conserve enthalpy

  IF (moist_has_happened == 1) THEN
     CALL conserve(nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, &
          p_full, p_half, p_surf, theta_final)
  ELSE
     CALL conserve(nlevels, theta_start, theta_adj_rlx, theta_adj_rlx, p_full, &
          p_half, p_surf, theta_final)
  END IF
  
  

END IF

! Note that if no convection happens, then theta_final, q_final and
! rainsec are left as they are initialised in 'initialise', so in this
! case theta_final = theta_start, q_final = q_start and rainsec = 0.

END SUBROUTINE convection_scheme

! ----------------------------------------------------------------------
! -- END OF CONVECTION SCHEME SUBROUTINE -- !
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! -- ROUTINES CALLED BY THE CONVECTION SCHEME -- !
! ----------------------------------------------------------------------

SUBROUTINE initialise(theta_start, theta_adjust,                        &
theta_noq_adjust, theta_adj_rlx, theta_noq_adj_rlx, theta_final,        &
q_start, q_adjust, q_adj_rlx, q_conserve, q_dummy, q_final,             &
rain_sec, final_level, moist_has_happened,                              &
convection_flag)

! 1) INITIALISE VARIABLES
! This subroutine initialises working arrays by setting them equal to
! input arrays.
!
! Further, the convbase and cloudbase flags are set to -1, as no
! convection has happened yet.
!
! Finally, rainsec is set to 0.

IMPLICIT NONE

! Inputs
REAL, INTENT(IN) :: theta_start(:), q_start(:)

! Outputs
REAL, INTENT(OUT) :: theta_adjust(:), theta_noq_adjust(:),              &
     theta_adj_rlx(:), theta_noq_adj_rlx(:), theta_final(:)
REAL, INTENT(OUT) :: q_adjust(:), q_adj_rlx(:), q_conserve(:),          &
     q_dummy(:), q_final(:)
REAL, INTENT(OUT) :: rain_sec
INTEGER, INTENT(OUT) :: final_level, moist_has_happened, convection_flag

theta_adjust = theta_start
theta_noq_adjust = theta_start
theta_adj_rlx = theta_start
theta_noq_adj_rlx = theta_start
theta_final = theta_start

q_adjust = q_start
q_adj_rlx = q_start
q_conserve = q_start
q_final = q_start
q_dummy = q_start

rain_sec = 0.0

final_level = 1
moist_has_happened = -1
convection_flag = -1

END SUBROUTINE initialise

! ----------------------------------------------------------------------

SUBROUTINE diagnose(nlevels, current_level, theta_start, q_start,       &
     p_full, p_surf, convbase, cloudbase, final_level, convection_flag, &
     moist_trigger, rh_use)

! 2) DIAGNOSE CONVECTION
! This subroutine decides whether it is possible for convection to
! originate on the model level under consideration, and if so, whether
! at any point above this level moist convection is possible (I.e
! whether a parcel rising from this level would saturate)
!
! The criterion for whether mixing is possible is the existence of a dry
! instability, that is, where d_theta / d_z is negative.
!
! Further the requirement for a rising parcel to saturate is
! q(convbase) > rh_crit * q_sat(level where it saturates).
!
! Note that rh_crit is a parameter of the scheme defined at the start.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels, current_level
REAL, INTENT(IN) :: theta_start(:), q_start(:), p_full(:)
REAL, INTENT(IN) :: p_surf

! Outputs
INTEGER, INTENT(INOUT) :: convection_flag
INTEGER, INTENT(OUT) :: convbase, cloudbase, final_level, moist_trigger
REAL, INTENT(OUT) :: rh_use

! For calculation
INTEGER :: j
REAL :: temperature, e_sat, q_sat, q_crit, rh

! initialise cloudbase
cloudbase = -1
convbase = -1
final_level = -1
moist_trigger = -1
rh_use = rh_crit

temperature = theta_start(current_level) * (p_full(current_level)       &
/p_zero) ** (kappa)
! DEPENDS ON: qsat_wat_mix
CALL qsat_wat_mix(q_sat,temperature,p_full(current_level),1,.FALSE.)
q_crit = rh_crit * q_sat
IF ((theta_start(current_level)  >=  theta_start(current_level+1)) .OR. &
(q_start(current_level)  >=  q_crit)) THEN
  
  IF (q_start(current_level)  >=  q_crit) THEN
     moist_trigger = 1
  END IF

  IF (q_start(current_level)  >=  0.0) THEN

    convbase = current_level

    IF (convection_flag  ==  -1) THEN ! Convection happens, change flag.
      convection_flag = 1
    END IF

    DO j = current_level , nlevels
      IF (j  ==  nlevels) THEN
        EXIT ! I.e no moist, so exit.
      END IF
      temperature = theta_start(j) * (p_full(j)/p_zero)**(Kappa)
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(q_sat,temperature,p_full(j),1,.FALSE.)
      q_crit = rh_crit * q_sat
      IF (q_start(j)  >=  q_crit) THEN
        cloudbase = j
        IF (cloudbase  ==  convbase) THEN
          cloudbase = convbase + 1 ! convbase cannot be cloudbase
        END IF
        EXIT
      END IF
   END DO

   ! Allows reference profile RH to adopt RH of layer where moist
   ! convection begins.
   
   IF (cloudbase /= -1) THEN
      temperature = theta_start(cloudbase) * (p_full(cloudbase)/p_zero)**(Kappa)
      CALL qsat_wat_mix(q_sat,temperature,p_full(cloudbase),1,.FALSE.)
      rh = q_start(cloudbase) / q_sat
      IF ((rh > rh_use).and.(rh<=1)) then
         rh_use = rh
      END IF
      IF (rh > 1) THEN
         rh_use = 1
      END IF
   END IF
   

  ELSE
    final_level = current_level + 1
  END IF
ELSE
  final_level = current_level +1
END IF

END SUBROUTINE diagnose

! ----------------------------------------------------------------------

SUBROUTINE convection(nlevels, theta_adjust, theta_noq_adjust, q_adjust,&
q_dummy, p_full, p_surf, convbase, cloudbase, final_level,              &
moist_has_happened, moist_trigger, rh_use)

! 3) DO CONVECTION
! This routine takes the values assigned to convbase and cloudbase and
! calls the appropriate routines for dry or moist convection.
!
! Once a convective event has finished it will record the level where it
! finished as final_level, which will then in turn be used by the
! routine 'diagnose' to find the values for convbase and cloudbase for
! any subsequent convective event.
!
! For the purpose of moist_has_happened, see section 5) (moist_convection)

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels, convbase, cloudbase
REAL, INTENT(IN) :: p_full(:)
REAL, INTENT(IN) :: p_surf, rh_use

! Outputs
REAL, INTENT(INOUT) :: theta_adjust(:), theta_noq_adjust(:),            &
     q_adjust(:), q_dummy(:)
INTEGER, INTENT(INOUT) :: moist_has_happened, moist_trigger
INTEGER, INTENT(OUT) :: final_level


IF (cloudbase  ==  -1) THEN
  CALL dry_convection(nlevels, convbase, theta_adjust,                  &
       theta_noq_adjust, q_adjust, final_level)
ELSE
  CALL moist_convection(nlevels, theta_adjust, theta_noq_adjust,        &
  q_adjust, q_dummy, p_full, p_surf, convbase, cloudbase, final_level,  &
  moist_has_happened, moist_trigger, rh_use)
END IF



END SUBROUTINE convection

! ----------------------------------------------------------------------

SUBROUTINE dry_convection(nlevels, convbase, theta_adjust,              &
theta_noq_adjust, q_adjust, final_level)

! 4) DRY CONVECTION
! This routine increases the potential temperature on model levels to
! those given by the dry adiabat, whilst the presence of a dry instability
! (as defined before) permits it.
!
! We also perform the same adjustment for theta_noq_adjust, a dummy 'dry
! convection only' profile that we use for enthalpy conservation.
!
! Specific humidity is also mixed in a similar way, by setting
! q_adjust(level_under_consideration) to q_start(convbase).
!
! Both of these adjustments will then be reduced to account for
! timestep / timescale.
!
! We also output 'final_level', the level where the convective event ends.
! It is from this level that the routine 'diagnose' will be re-called, in
! an effort to find additional thermal instabilities that would permit
! convection.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels, convbase

! Outputs
REAL, INTENT(INOUT) :: theta_adjust(:), theta_noq_adjust(:),            &
     q_adjust(:)
INTEGER, INTENT(OUT) :: final_level

! Calculation
INTEGER :: j

final_level = -1 ! if this remains unchanged then, in circumstances
                 ! where dry_convection has been called by moist_
                 ! convection, it signals to moist_convection that
                 ! the event is able to proceed. I.e the dry mixing
                 ! reached the level cloudbase.

DO j = convbase+1 , nlevels
  IF ((theta_adjust(j-1)  >=  theta_adjust(j)) .AND.                    &
  (q_adjust(j-1)  >=  0.0)) THEN
    theta_adjust(j) = theta_adjust(j-1)
    theta_noq_adjust(j) = theta_noq_adjust(j-1)
    IF ((q_adjust(j-1)  >=  0.0) .AND. (q_adjust(j)  <                  &
    q_adjust(j-1))) THEN
      q_adjust(j) = q_adjust(j-1)
    END IF
  ELSE
    final_level = j
    EXIT
  END IF
END DO

END SUBROUTINE dry_convection

! ----------------------------------------------------------------------

SUBROUTINE moist_convection(nlevels, theta_adjust, theta_noq_adjust,    &
q_adjust, q_dummy, p_full, p_surf, convbase, cloudbase, final_level,    &
moist_has_happened, moist_trigger, rh_use)

! 5) MOIST CONVECTION
! This routine is called if we know that a level cloudbase exists
! such that a parcel of air rising from convbase would saturate.
!
! In this situation the theta profile is adjustted beyond that of a
! dry reference profile to the reference profile for the moist
! pseudo-adiabat. This is as a result of the latent heating provided
! by saturation.
!
! The specific humidity within the columns is mixed as before, unless
! this would lead to a level having specific humidity greater than
! saturation specific humidity, in which case we set specific humidity
! equal to saturation specific humidity.
!
! If cloudbase-convbase >= 2 then we perform dry convection until
! we reach the level cloudbase.
!
! If dry convection ceases to be possible before we reach the level
! cloudbase then we claim that a parcel rising from convbase does not
! reach a level at which it saturates, and so we end the convective event
! accordingly without actually performing moist convection. Moist
! covection may still be possible in such a scenario if, when we re-call
! diagnose from 'final_level', it finds another level convbase with
! accompanying level cloudbase.
!
! Finally it returns an integer, moist_has_happened, which signals whether
! moist convection has infact taken place. This is initalised in the
! routine 'initialise' to be moist_has_happened = -1. If at any point
! moist convection takes place this is changed to 1.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels, convbase, cloudbase
REAL, INTENT(IN) :: p_full(:)
REAL, INTENT(IN) :: p_surf, rh_use

! Outputs
INTEGER, INTENT(OUT) :: final_level
INTEGER, INTENT(INOUT) :: moist_has_happened, moist_trigger
REAL, INTENT(INOUT) :: theta_adjust(:), theta_noq_adjust(:),            &
     q_adjust(:), q_dummy(:)

! For calcualtion
INTEGER :: j
REAL :: temperature(nlevels)
REAL :: e_sat, q_sat, r_sat, virtual_temperature, gamma_s, gamma_p, rho, &
     q_adjust_save

final_level = -1 ! final_level = -1 allows moist convection to happen.

CALL calc_temp(theta_adjust, p_full, temperature)

IF ((cloudbase-convbase)  >=  2) THEN
  ! Note that in this situation we only want to perform dry convection
  ! up until the level below cloudbase, so we use cloudbase-1 as the
  ! 'nlevels' input.
  CALL dry_convection((cloudbase-1), convbase, theta_adjust,            &
  theta_noq_adjust, q_adjust, final_level)
  CALL calc_temp(theta_adjust, p_full, temperature)
END IF

                              ! Dry convection was able to reach level
IF (final_level  ==  -1) THEN ! below cloudbase, so we now do moist
                              ! convection.
  DO j = cloudbase , nlevels
    IF ((theta_adjust(j-1)  >=  theta_adjust(j)) .OR.                   &
         (moist_trigger  ==  1)) THEN

      IF (moist_trigger  ==  1) THEN
        moist_trigger = -1 ! reset moist trigger
      END IF

      IF (q_adjust(j-1)  >=  0.0) THEN

        ! DEPENDS ON: qsat_wat_mix
         CALL qsat_wat_mix(q_sat,temperature(j),p_full(j),1,.FALSE.)
         q_sat = rh_use * q_sat

        IF (q_adjust(j)  <   q_sat) THEN
           q_adjust_save = q_adjust(j)
          q_adjust(j) = q_sat ! adjust SH profile
        END IF

        ! Now we adjust potential temperature
        ! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(r_sat,temperature(j),p_full(j),1,.TRUE.)
        r_sat = rh_use * r_sat
        gamma_s = (-g) * ((1+r_sat)*(1 + (lc*r_sat)/                    &
        (r*temperature(j-1))))   /   (cp+r_sat*hcapv+(repsilon+r_sat)*  &
        ((r_sat*lc**2)/(r*temperature(j-1)**2))) ! sat dT/dz
        ! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(r_sat,temperature(j-1),p_full(j-1),1,.TRUE.)
        r_sat = rh_use * r_sat
        virtual_temperature = temperature(j-1) * (1+r_sat/repsilon) / (1+r_sat)
        ! virtual_temperature used for density of moist air
        rho = p_full(j-1) / (r*virtual_temperature) ! moist air density
        gamma_p = (1.0/(rho*g)) * gamma_s ! saturated dT/dp
        IF ((temperature(j-1)+(gamma_p*(p_full(j-1)-p_full(j)))) <= &
             temperature(j)) THEN
           final_level = j
           q_adjust(j) = q_adjust_save
           EXIT
        END IF
        
 
        temperature(j) = temperature(j-1) + (gamma_p *                  &
        (p_full(j-1)-p_full(j))) ! this is reference profile temperature
        theta_adjust(j) = temperature(j) * ((p_zero/p_full(j)) **       &
        (kappa)) ! adjust potential temperature accordingly

        IF (moist_has_happened  ==  -1) THEN ! set moist_has_happened=1
          moist_has_happened = 1
        END IF

        ! Adjust dummy dry potential temperature profile
        IF (theta_noq_adjust(j-1)  >=  theta_noq_adjust(j)) THEN
          theta_noq_adjust(j) = theta_noq_adjust(j-1)
        END IF

      ELSE
        final_level = j
        EXIT
      END IF

    ELSE ! record final level
      final_level = j
      EXIT
    END IF
  END DO
END IF

END SUBROUTINE moist_convection

! ----------------------------------------------------------------------

SUBROUTINE relax(theta_start, theta_adjust, theta_adj_rlx, q_start,     &
q_adjust, q_adj_rlx)

! 6) TIMESCALE
! This routine modifies the amount by which potential temperature and
! specific humidity are increased to account for convective mixing
! timescale and model timestep.
!
! In cases where moist convection has happened, this routine will be
! called twice; once for theta_adjust, and once for theta_noq_adjust.

IMPLICIT NONE

! Input
REAL, INTENT(IN) :: theta_start(:), theta_adjust(:), q_start(:),        &
     q_adjust(:)

! Output
REAL, INTENT(INOUT) :: theta_adj_rlx(:), q_adj_rlx(:)

theta_adj_rlx = theta_start + ((theta_adjust-theta_start)) *            &
(1.0 - EXP(-timestep/timescale))
q_adj_rlx = q_start + ((q_adjust-q_start)) * (1.0 - EXP(-timestep/      &
timescale))

END SUBROUTINE relax

! ----------------------------------------------------------------------

SUBROUTINE moisture_conserve(nlevels, q_start, q_adj_rlx, p_full,       &
p_half, p_surf, q_conserve)

! 7) MOISTURE CONVSERVE
! This routine ensures that at this point in the scheme there is no net
! change in moisture within the column (this will change later, as some
! is lost through rainfall).
!
! To do this, the routine find_bounds is used to find the top and
! bottom of convective events that have occured within the column.
!
! For each event, we find the total mass of water contained on the levels
! involved in this event both before and after convection. We then take
! the ratio of water : before/after and multiply the specfic humidity on
! each level involved in the event by this ratio.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
REAL, INTENT(IN) :: q_start(:), q_adj_rlx(:), p_full(:)
REAL, INTENT(IN) :: p_half(:)
REAL, INTENT(IN) :: p_surf

! Outputs
REAL, INTENT(OUT) :: q_conserve(:)

! For calculation
INTEGER :: stopnow, bottom, top, j, k ,n ! stopnow, bottom, top and are
                                         ! used by find_bounds
REAL :: p_change(nlevels), q_start_mass(nlevels),                       &
     q_adj_rlx_mass(nlevels), q_inc(nlevels)
REAL :: q_start_total_mass, q_adj_rlx_total_mass, q_inc_ratio

stopnow = 0 ! Initialise stopnow
q_conserve = q_adj_rlx

! Find p_change; change in pressure from bottom to top of levels
p_change = p_half(1:nlevels) - p_half(2:(nlevels+1))

! Find increase in specific humidity after convection
! find_bounds uses the difference between the start and adjusted profile
! to determine where a convective event begins and ends. For logistical
! reasons it requires a positive difference. For more information, see
! description of find_bounds.
q_inc = q_adj_rlx - q_start
DO j = 1 , nlevels
  q_inc(j) = ABS(q_inc(j))
END DO

n = 2 ! part of find_bounds set up.

DO j = 1, nlevels

  IF (j  ==  n) THEN
    CALL find_bounds(n, bottom, top, q_inc, stopnow, nlevels)
    IF (stopnow  /=  1) THEN

      DO k = bottom , top ! find masses on each level within event
        q_start_mass(k) = 0.0
        q_adj_rlx_mass(k) = 0.0
        q_start_mass(k) = ABS(q_start(k)) * p_change(k) / g
        q_adj_rlx_mass(k) = ABS(q_adj_rlx(k)) * p_change(k) / g
      END DO

      ! find total masses for each event
      q_start_total_mass = SUM(q_start_mass(bottom:top))
      q_adj_rlx_total_mass = SUM(q_adj_rlx_mass(bottom:top))

      ! find ratio before / after
      q_inc_ratio = q_start_total_mass / q_adj_rlx_total_mass

      ! set q_conserve accordingly
      DO k = bottom , top
        q_conserve(k) = ABS(q_adj_rlx(k)) * q_inc_ratio
      END DO

    ELSE
      EXIT
    END IF
  END IF

END DO

END SUBROUTINE moisture_conserve

! ----------------------------------------------------------------------

SUBROUTINE rain(nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, &
     q_conserve, p_full, p_half, p_surf, q_final, rain_sec)

! 8) RAIN / LATENT
! This routine calculates the total amount of latent heating within a
! convective event that was required for the moist convection to take
! place.
!
! It then removes the mass of water necessary to do this in the form of
! instantaneous precipitation and adjusts the specific humidity of the
! levels under consideration accordingly.
!
! Total rainfall in the column is found via Q = mass * L where Q is
! latent heating.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
REAL, INTENT(IN) :: theta_start(:), theta_noq_adj_rlx(:),               &
     q_conserve(:), p_full(:), p_half(:)
REAL, INTENT(IN) :: p_surf

! Outputs
REAL, INTENT(INOUT) :: theta_adj_rlx(:)
REAL, INTENT(OUT) :: q_final(:)
REAL, INTENT(OUT) :: rain_sec

! For calculation
INTEGER :: j, k, n, stopnow, bottom, top
REAL :: q_level(nlevels), latent(nlevels), temperature(nlevels),        &
     p_change(nlevels), temperature_noq(nlevels), theta_inc(nlevels)
REAL :: event_q_total, total_latent, total_latent_new, rain_out,        &
rain_out_new, latent_ratio, temperature_new

stopnow  = 0 ! initialise stop now
rain_sec = 0.0

! Find p_change; change in pressure from bottom to top of levels
p_change = p_half(1:nlevels) - p_half(2:(nlevels+1))

q_final = q_conserve

! Find increase in potential temperature after convection
! find_bounds uses the difference between the start and adjusted profile
! (as before.. see description in moisture_conserve)
theta_inc = theta_adj_rlx - theta_start
DO j = 1 , nlevels
  theta_inc(j) = ABS(theta_inc(j))
END DO

n = 2 ! part of find_bounds set up.

! calculate temperature array for dummy 'only dry convection' theta
! profile
CALL calc_temp(theta_noq_adj_rlx, p_full, temperature_noq)

DO j = 1 , nlevels

  IF (j  ==  n) THEN
    CALL find_bounds(n, bottom, top, theta_inc, stopnow, nlevels)

    IF (stopnow  /=  1) THEN

      ! calculate temperature array for actual theta profile
      CALL calc_temp(theta_adj_rlx, p_full, temperature)

      ! reset all arrays and totals to 0 for new event
      DO k = 1 , nlevels
        latent(k) = 0.0
        q_level(k) = 0.0
      END DO
      event_q_total = 0.0
      total_latent = 0.0
      total_latent_new = 0.0
      latent_ratio = 0.0
      rain_out = 0.0
      rain_out_new = 0.0

      ! Rainfall / latent heating adjustment

      DO k = bottom , top
        latent(k) = p_change(k)/g * cp * (temperature(k) -              &
        temperature_noq(k))
      END DO

      ! Calculate rainfall within this column from the amount
      ! of latent heating that has taken place.
      total_latent  = SUM(latent(bottom:top))
      IF (total_latent  >   0.0) THEN
        rain_out = total_latent / lc
        rain_sec = rain_sec + (rain_out / timestep )

        ! find total water on levels involved in convective event.
        DO k = bottom , top
          q_level(k) = (p_change(k) / g) * q_conserve(k)
        END DO
        event_q_total = SUM(q_level(bottom:top))

        ! The following code handles situations where the amount of rain
        ! demanded by a convective event is greater than the amount of
        ! rain stored in the levels participating in the convective event.
        ! In this scenario we limit rain_out to the amount of moisture
        ! present in the event and revise our convective potential
        ! temperature adjustment accordingly.
        IF (rain_out  >   event_q_total) THEN
          rain_out_new = event_q_total
          rain_sec = rain_sec - (rain_out/timestep) + (rain_out_new/    &
          timestep)
          total_latent_new = rain_out_new * lc
          latent_ratio  = total_latent_new / total_latent
          latent = latent * latent_ratio
          DO k = bottom , top
            temperature_new = temperature_noq(k) + latent(k)/           &
            ((p_change(k)/g)*cp)
            theta_adj_rlx(k) = temperature_new / ((p_full(k)/           &
            p_zero)**(kappa))
          END DO
          event_q_total = rain_out_new
          rain_out = rain_out_new
        END IF

         ! Resume rainfall / latent heating adjustment

         ! Rescale specific humidity on levels within this event to
         ! account for rainfall
        DO k = bottom , top
          q_level(k) = q_level(k) * (event_q_total - rain_out) /        &
          event_q_total
          q_final(k) = (g/p_change(k)) * q_level(k)
        END DO
      END IF
    ELSE
      EXIT
    END IF
  END IF

END DO

END SUBROUTINE rain

! ----------------------------------------------------------------------

SUBROUTINE conserve(nlevels, theta_start, theta_adj_rlx,                &
theta_noq_adj_rlx, p_full, p_half, p_surf, theta_final)

! 9) ENTHALPHY CONSERVE
! This routine conserves energy that has been added to the column in the
! form of dry heating.
!
! Within a given convective event between the levels bottom and top
! (found by find_bounds), we find out how much energy was required to
! allow the increase in potential temperature from theta_start to
! theta_noq_adj_rlx (our dummy dry profile).
!
! In order to calculate how much energy we have added during the event
! under consideration we use dQ = cp * dT * mass.

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
REAL, INTENT(IN) :: theta_start(:), theta_adj_rlx(:),                   &
     theta_noq_adj_rlx(:), p_full(:)
REAL, INTENT(IN) :: p_half(:)
REAL, INTENT(IN) :: p_surf

! Output
REAL, INTENT(OUT) :: theta_final(:)

! For calculation
REAL :: temperature_start(nlevels), temperature_adj_rlx(nlevels),       &
     temperature_noq_adj_rlx(nlevels), temperature_final(nlevels),      &
     dry_heating(nlevels), mass(nlevels), p_change(nlevels),            &
     theta_dry_inc(nlevels)
REAL :: total_heating, total_mass, temperature_correction
INTEGER :: stopnow, bottom, top, n, j, k

! Find temperatures
CALL calc_temp(theta_start, p_full, temperature_start)
CALL calc_temp(theta_adj_rlx, p_full, temperature_adj_rlx)
CALL calc_temp(theta_noq_adj_rlx, p_full, temperature_noq_adj_rlx)
CALL calc_temp(theta_final, p_full, temperature_final)

stopnow = 0 ! initialise stopnow

! Find p_change
p_change = p_half(1:nlevels) - p_half(2:(nlevels+1))

theta_dry_inc = theta_noq_adj_rlx - theta_start
theta_final = theta_adj_rlx

n = 2 ! required for find_bounds

DO j = 1 , nlevels

  IF (j  ==  n) THEN
    CALL find_bounds(n, bottom, top, theta_dry_inc, stopnow, nlevels)
    IF (stopnow  /=  1) THEN
      DO k = bottom , top
        mass(k) = 0.0
        dry_heating(k) = 0.0
        mass(k) = p_change(k) / g
        dry_heating(k) = cp * (temperature_noq_adj_rlx(k)               &
        - temperature_start(k)) * mass(k)
      END DO
      total_mass = SUM(mass(bottom:top))
      total_heating = SUM(dry_heating(bottom:top))
      temperature_correction = total_heating / (cp * total_mass)
      DO k = bottom , top
        temperature_final(k) = temperature_adj_rlx(k)                   &
        - temperature_correction
        theta_final(k) = temperature_final(k) / ((p_full(k)/p_zero)     &
        ** (kappa))
      END DO
    ELSE
      EXIT
    END IF
  END IF
END DO

END SUBROUTINE conserve

! ----------------------------------------------------------------------
! -- END OF ROUTINES CALLED BY CONVECTION SCHEME -- !
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! -- EXTRA ROUTINES -- !
! ----------------------------------------------------------------------

SUBROUTINE find_bounds(n, bottom, top, inc, stopnow, nlevels)

 ! FIND BOUNDS
 ! This routine finds the bounds for a convective event from either the
 ! potential temperature or specific humidity increment for that event.
 !
 ! Note that this routine also outputs the integer n which is the level
 ! above the top level for the convective event under consideration. In
 ! other words, this is the first level above a convective event where a
 ! different convective event could be found.

IMPLICIT NONE

REAL, INTENT(IN) :: inc(:)
INTEGER, INTENT(IN) :: nlevels
INTEGER :: j
INTEGER, INTENT(INOUT) :: stopnow, n
INTEGER, INTENT(OUT) :: bottom, top

bottom = -1

DO j = n,nlevels

  ! find bottom
  IF ((bottom  ==  -1) .AND. (inc(j)  >   0)) THEN
    bottom = j-1
  ELSE IF ((bottom  ==  -1) .AND. (j  ==  nlevels)) THEN
    stopnow = 1
    EXIT
  END IF

  ! find top
  IF (j  ==  nlevels) THEN
    top = j-1
    EXIT
  ELSE IF (((j-1)  /=  bottom) .AND. (inc(j)  ==  0)) THEN
    IF (bottom  /=  -1) THEN
      top = j-1
      EXIT
    END IF
  END IF

END DO

n = top + 1

END SUBROUTINE find_bounds

! ----------------------------------------------------------------------

SUBROUTINE calc_temp(theta, p_full, temp)

! CALCULATE TEMPERATURE
! This routine calculates and array containing temperatures from an
! array containing potential temperatures.

IMPLICIT NONE

REAL, INTENT(IN) :: theta(:), p_full(:)
REAL, INTENT(OUT) :: temp(:)

temp = theta * (p_full/p_zero)**(kappa)

END SUBROUTINE calc_temp

! ----------------------------------------------------------------------
! -- END OF EXTRA ROUTINES -- !
! ----------------------------------------------------------------------

!-------    
!-- QSAT_WAT_MIX FROM UM
!-------

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Saturation Specific Humidity Scheme (Qsat_Wat): Vapour to Liquid.
! Subroutine Interface:
SUBROUTINE qsat_wat_mix (                                         &
!      Output field
        QmixS                                                           &
!      Input fields
      , t, p                                                            &
!      Array dimensions
      , npnts                                                           &
!      logical control
      , lq_mix                                                          &
        )


!USE vectlib_mod, ONLY: oneover_v
!USE yomhook, ONLY: lhook, dr_hook
!USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Purpose:
!   Returns a saturation specific humidity or mixing ratio given a
!   temperature and pressure using the saturation vapour pressure
!   calculated using the Goff-Gratch formulae, adopted by the WMO as
!   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
!   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
!   Phyiscal and Chemical properties or air, P35.
!
!   Values in the lookup table are over water above and below 0 deg C.
!
!   Note : For vapour pressure over water this formula is valid for
!   temperatures between 373K and 223K. The values for saturated vapour
!   over water in the lookup table below are out of the lower end of
!   this range. However it is standard WMO practice to use the formula
!   below its accepted range for use with the calculation of dew points
!   in the upper atmosphere
!
! Method:
!   Uses lookup tables to find eSAT, calculates qSAT directly from that.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service
!
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29
!
! Declarations:
!
!  Global Variables:----------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
!
! arguments with intent in. ie: input variables.

INTEGER, INTENT(IN) :: npnts
! Points (=horizontal dimensions) being processed by qSAT scheme.

REAL, INTENT(IN)  :: t !  Temperature (K).
REAL, INTENT(IN)  :: p !  Pressure (Pa).

LOGICAL, INTENT(IN)  :: lq_mix
              !  .true. return qsat as a mixing ratio
              !  .false. return qsat as a specific humidity

! arguments with intent out

REAL, INTENT(OUT)   ::  QmixS
       ! Output Saturation mixing ratio or saturation specific
       ! humidity at temperature T and pressure P (kg/kg).

!-----------------------------------------------------------------------
!  Local scalars
!-----------------------------------------------------------------------

INTEGER :: itable    ! Work variable

REAL :: atable       ! Work variable

REAL :: fsubw
      ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A PURE
      ! WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR.

REAL :: tt

REAL :: vect_in, vect_out
      ! temp array for calculating reciprocals

REAL, PARAMETER :: R_delta_T = 1.0/delta_T

INTEGER :: i

!INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!REAL(KIND=jprb)               :: zhook_handle

!CHARACTER(LEN=*), PARAMETER :: RoutineName='QSAT_WAT_MIX'

!-----------------------------------------------------------------------

!IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


! loop over points
DO i=1,npnts

  tt = MAX(T_low,t)
  tt = MIN(T_high,tt)
  atable = (tt - T_low + delta_T) * R_delta_T
  itable = atable
  atable = atable - itable

    !      Use the lookup table to find saturated vapour pressure, and store
    !      it in qmixs.

  QmixS   = (1.0 - atable)    * es(itable)    +              &
       atable*es(itable+1)

  IF (lq_mix) THEN

      !      Compute the factor that converts from sat vapour pressure in a
      !      pure water system to sat vapour pressure in air, fsubw.
      !      This formula is taken from equation A4.7 of Adrian Gill's book:
      !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
      !      of pressure in mb and temperature in Celsius, so conversion of
      !      units leads to the slightly different equation used here.

    fsubw    = 1.0 + 1.0e-8*p   * ( 4.5 +                    &
           6.0e-4*( t   - zerodegC ) * ( t - zerodegC ) )

      !      Multiply by fsubw to convert to saturated vapour pressure in air
      !      (equation A4.6 of Adrian Gill's book).

    QmixS   = QmixS   * fsubw

      !      Now form the accurate expression for qmixs, which is a rearranged
      !      version of equation A4.3 of Gill's book.

      !-----------------------------------------------------------------------
      ! For mixing ratio,  rsat = epsilon *e/(p-e)
      ! e - saturation vapour pressure
      ! Note applying the fix to qsat for specific humidity at low pressures
      ! is not possible, this implies mixing ratio qsat tends to infinity.
      ! If the pressure is very low then the mixing ratio value may become
      ! very large.
      !-----------------------------------------------------------------------

    vect_in = ( MAX(p,  1.1*QmixS)   - QmixS )

  ELSE

      !      Compute the factor that converts from sat vapour pressure in a
      !      pure water system to sat vapour pressure in air, fsubw.
      !      This formula is taken from equation A4.7 of Adrian Gill's book:
      !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
      !      of pressure in mb and temperature in Celsius, so conversion of
      !      units leads to the slightly different equation used here.

    fsubw    = 1.0 + 1.0e-8*p   * ( 4.5 +                    &
           6.0e-4*( t   - zerodegC ) * ( t - zerodegC ) )

      !      Multiply by fsubw to convert to saturated vapour pressure in air
      !      (equation A4.6 of Adrian Gill's book).

    QmixS   = QmixS   * fsubw

      !      Now form the accurate expression for qmixs, which is a rearranged
      !      version of equation A4.3 of Gill's book.

      !-----------------------------------------------------------------------
      ! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
      !
      ! Note that at very low pressure we apply a fix, to prevent a
      ! singularity (qsat tends to 1. kg/kg).
      !-----------------------------------------------------------------------
    vect_in = (MAX(p, QmixS)-one_minus_epsilon*QmixS)

  END IF
END DO

DO i = 1 , npnts
   vect_out = 1 / vect_in
END DO

DO i=1, npnts
  QMixS = repsilon*QmixS * vect_out
END DO

!IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
!RETURN
END SUBROUTINE qsat_wat_mix

!------
!-- END OF QSAT_WAT_MIX
!--------

END MODULE llcs
