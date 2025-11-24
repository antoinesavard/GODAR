subroutine forcing (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_forcings.h"
    include "CB_diagnostics.h"
    include "CB_mpi.h"

    integer, intent(in) :: i

    double precision :: fdax, fday, fsax, fsay
    double precision :: fdwx, fdwy, fswx, fswy

    double precision :: log_profile, L2norm
    double precision :: shelter_coeff_a, shelter_coeff_w

    double precision :: atm_vel_norm, ocn_vel_norm, atm_log, ocn_log

    ! unitless minimum sheltering coefficient
    shelter_coeff_a = hsfa_min_r(i)
    shelter_coeff_w = hsfw_min_r(i)

    ! local variables to speed up the code
    atm_vel_norm = L2norm(ua - u(i), va - v(i))
    ocn_vel_norm = L2norm(uw - u(i), vw - v(i))
    atm_log = log_profile(hfa(i), z0w)
    ocn_log = log_profile(hfw(i), z0w)

    !###############################################################
    if ( atm_vel_norm .eq. 0d0 ) then
        fdax = 0d0
        fday = 0d0
        fsax = 0d0
        fsay = 0d0
    else
        ! wind drag forcing
        fdax     = rhoair * Cdair * hfa(i) * r(i) * (ua - u(i)) *   &
                    atm_vel_norm * atm_log * shelter_coeff_a *   &
                    pi / 2d0

        fday     = rhoair * Cdair * hfa(i) * r(i) * (va - v(i)) *   &
                    atm_vel_norm * atm_log * shelter_coeff_a *   &
                    pi / 2d0

        ! wind skin forcing
        fsax     = rhoair * pi * r(i) ** 2 * Csair *    &
                    (ua - u(i)) / atm_vel_norm * (      &
                    atm_vel_norm ** 2 + r(i) ** 2       &
                    * ABS(omega(i)) ** 2 / 4d0)

        fsay     = rhoair * pi * r(i) ** 2 * Csair *    &
                    (va - v(i)) / atm_vel_norm * (      &
                    atm_vel_norm ** 2 + r(i) ** 2       &
                    * ABS(omega(i)) ** 2 / 4d0)
    end if

    !###############################################################
    if ( ocn_vel_norm .eq. 0d0 ) then
        fdwx = 0d0
        fdwy = 0d0
        fswx = 0d0
        fswy = 0d0
    else
        ! water drag forcing
        fdwx = rhowater * Cdwater * hfw(i) * r(i) * (uw - u(i)) *   &
                ocn_vel_norm * ocn_log * shelter_coeff_w *       &
                pi / 2d0

        fdwy = rhowater * Cdwater * hfw(i) * r(i) * (vw - v(i)) *   &
                ocn_vel_norm * ocn_log * shelter_coeff_w *       &
                pi / 2d0

        ! water skin forcing
        fswx = rhowater * pi * r(i) ** 2 * Cswater *    &
                (uw - u(i)) / ocn_vel_norm * (          &
                ocn_vel_norm ** 2 + r(i) ** 2           &
                * ABS(omega(i)) ** 2 / 4d0)

        fswy = rhowater * pi * r(i) ** 2 * Cswater *    &
                (vw - v(i)) / ocn_vel_norm * (          &
                ocn_vel_norm ** 2 + r(i) ** 2           &
                * ABS(omega(i)) ** 2 / 4d0)
    end if

    !###############################################################
    ! total forcing from air and water
    fax(i) = fdax + fsax
    fay(i) = fday + fsay

    fwx(i) = fdwx + fswx
    fwy(i) = fdwy + fswy

    ! torque induced drag due to rotation of floes when no speed
	! if speed, use second expression valid for |U| >> |omega*r|
	if ( atm_vel_norm < 1d-2 ) then
		ma(i)  = - 2d0 * pi / 5d0 * r(i) ** 5 * rhoair * Csair * &
                omega(i) * ABS(omega(i))
	else
		ma(i) = - 3d0 / 4d0 * pi * rhoair * Csair * &
				sqrt( ua ** 2 + va ** 2) * omega(i) * r(i) ** 4
	end if

	if ( ocn_vel_norm < 1d-2 ) then
    	mw(i) = - 2d0 * pi / 5d0 * r(i) ** 5 * rhowater * Cswater * &
				omega(i) * ABS(omega(i))
	else
		mw(i) = - 3d0 / 4d0 * pi * rhowater * Cswater * &
				sqrt( uw ** 2 + vw ** 2) * omega(i) * r(i) ** 4
	end if

    ! compute the stress using cauchy stress formula due to the winds and currents
    ! off diag are always 0
    sigxx_aw(i) = (fax(i) + fwx(i)) * r(i)
    sigyy_aw(i) = (fay(i) + fwy(i)) * r(i)
    sigxy_aw(i) = 0d0
    sigyx_aw(i) = 0d0

end subroutine forcing


subroutine sheltering (j, i)

    implicit none

    include "parameter.h"
    include "CB_forcings.h"
    include "CB_variables.h"

    integer, intent(in) :: i, j

    double precision :: S_shelter

    ! sheltering height from air and water
    hsfa(i, j) = S_shelter(hfa(i), hfa(j), deltan(j,i), cosa(j,i), &
                    sina(j,i), ua, va)
    hsfw(i, j) = S_shelter(hfw(i), hfw(j), deltan(j,i), cosa(j,i), &
                    sina(j,i), uw, vw)

    ! sheltering for the reverse direction
    hsfa(j, i) = S_shelter(hfa(j), hfa(i), deltan(j,i), -cosa(j,i), &
                    -sina(j,i), ua, va)
    hsfw(j, i) = S_shelter(hfw(j), hfw(i), deltan(j,i), -cosa(j,i), &
                    -sina(j,i), uw, vw)

end subroutine sheltering


subroutine coriolis (i)

	implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

	integer, intent(in) :: i
    double precision :: omega_earth

    omega_earth = 7.2921d-5

	fcorx(i)  = mass(i) * 2 * omega_earth * sin(lat) * v(i)

    fcory(i)  = - mass(i) * 2 * omega_earth * sin(lat) * u(i)

end subroutine coriolis


double precision function log_profile (hf, z0)

    implicit none

    double precision, intent(in) :: hf, z0

    ! this is the true profile to use but it is unnecessarily complex
    ! we can use a nice little approximation that is valid within 1%
    ! because z0 is very small in compared to hf
    ! log_profile = ( log( hf / z0 ) ** 2.0d0 - 2.0d0 *           &
    !                 log( hf / z0 ) + 2.0d0 - 2.0d0 * z0 / hf) / &
    !                 log( 10 / z0 ) ** 2.0d0

    log_profile = ( log( hf / z0 ) / log( 10 / z0 ) ) ** 2.0d0

end function log_profile

double precision function heaviside (x)

    implicit none

    double precision, intent(in) :: x

    ! note how we use the -fno-sign-zero compiler option to prevent 
    ! signed zeros coming from here
    heaviside = sign(5d-1, x) + 0.5

end function heaviside

double precision function S_shelter (hfi, hfj, deltan, cosa, sina, uf, vf)

    implicit none

    include "CB_const.h"

    double precision, intent(in) :: hfi, hfj, deltan
    double precision, intent(in) :: uf, vf, cosa, sina
    
    double precision :: heaviside, L2norm
    
    ! D: is the distance factor
    ! G: is the upwind factor
    ! costheta: where theta is the angle between the wind vector 
    ! and the vector joining the center of the disks.
    ! wnorm: is the L2norm of the winds
    double precision :: D, G, costheta, wnorm

    wnorm = L2norm(uf, vf) + 1d-20
    costheta = (uf * cosa + vf * sina) / wnorm

    D = 1 -                                                 &
        ( 1 - heaviside(deltan) ) *                         &
        ( 1 - exp(- 0.18 * abs(deltan / costheta) / hfj) )  &
        -                                                   &
        heaviside(deltan) *                                 &
        heaviside(hfi - hfj) * (1 - hfj / hfi)
    
    G = heaviside(costheta) * costheta

    S_shelter = 1 - D * G

end function S_shelter


double precision function L2norm (u, v)

    implicit none
    
    double precision, intent(in) :: u, v

    L2norm = sqrt( u ** 2d0 + v ** 2d0 )

end function L2norm
