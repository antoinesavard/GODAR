subroutine forcing (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_forcings.h"

    integer, intent(in) :: i

    double precision :: fdax, fday, fsax, fsay
    double precision :: fdwx, fdwy, fswx, fswy

    double precision :: log_profile, L2norm
    double precision :: shelter_coeff_a(n), shelter_coeff_w(n)

    ! unitless maximum sheltering
    shelter_coeff_a = maxval(hsfa, dim=1)
    shelter_coeff_w = maxval(hsfw, dim=1) 
    
    ! wind drag forcing
    fdax     = rhoair * Cdair * hfa(i) * r(i) * (ua - u(i)) *      &
                L2norm(ua - u(i), va - v(i)) * log_profile(hfa(i), &
                max(z0w, shelter_coeff_a(i) * hfa(i))) *           &
                shelter_coeff_a(i) * pi / 2d0

    fday     = rhoair * Cdair * hfa(i) * r(i) * (va - v(i)) *      &
                L2norm(ua - u(i), va - v(i)) * log_profile(hfa(i), &
                max(z0w, shelter_coeff_a(i) * hfa(i))) *           &
                shelter_coeff_a(i) * pi / 2d0

    ! wind skin forcing
    fsax     = rhoair * pi * r(i) ** 2 * Csair *               &
                (ua - u(i)) / L2norm(ua - u(i), va - v(i)) * ( &
                L2norm(ua - u(i), va - v(i)) ** 2 + r(i) ** 2  &
                * ABS(omega(i)) ** 2 / 4d0)

    fsay     = rhoair * pi * r(i) ** 2 * Csair *               &
                (va - v(i)) / L2norm(ua - u(i), va - v(i)) * ( &
                L2norm(ua - u(i), va - v(i)) ** 2 + r(i) ** 2  &
                * ABS(omega(i)) ** 2 / 4d0)

    ! water drag forcing
    fdwx     = rhowater * Cdwater * hfw(i) * r(i) * (uw -            &
                u(i)) * L2norm(uw - u(i), vw - v(i)) *               &
                log_profile(hfw(i), max(z0w, shelter_coeff_w(i) *    &
                hfw(i))) * shelter_coeff_w(i) * pi / 2d0

    fdwy     = rhowater * Cdwater * hfw(i) * r(i) * (vw -            &
                v(i)) * L2norm(uw - u(i), vw - v(i)) *               &
                log_profile(hfw(i), max(z0w, shelter_coeff_w(i) *    &
                hfw(i))) * shelter_coeff_w(i) * pi / 2d0

    ! water skin forcing
    fswx     = rhowater * pi * r(i) ** 2 * Cswater *               &
                (uw - u(i)) / L2norm(uw - u(i), vw - v(i)) * (     &
                L2norm(uw - u(i), vw - v(i)) ** 2 + r(i) ** 2      &
                * ABS(omega(i)) ** 2 / 4d0)

    fswy     = rhowater * pi * r(i) ** 2 * Cswater * &
                (vw - v(i)) / L2norm(uw - u(i), vw - v(i)) * (     &
                L2norm(uw - u(i), vw - v(i)) ** 2 + r(i) ** 2      &
                * ABS(omega(i)) ** 2 / 4d0)

    ! total forcing from air and water
    fax(i) = fdax + fsax
    fay(i) = fday + fsay

    fwx(i) = fdwx + fswx
    fwy(i) = fdwy + fswy

    ! torque induced drag due to rotation of floes when no speed
	! if speed, use second expression valid for |U| >> |omega*r|
	if ( L2norm(ua - u(i), va - v(i)) .eq. 0d0 ) then
        fax(i) = 0d0
        fay(i) = 0d0
		ma(i)  = - 2d0 * pi / 5d0 * r(i) ** 5 * rhoair * Csair * &
                omega(i) * ABS(omega(i))
	else
		ma(i) = - 3d0 / 4d0 * pi * rhoair * Csair * &
				sqrt( ua ** 2 + va ** 2) * omega(i) * r(i) ** 4
	end if

	if ( L2norm(uw - u(i), vw - v(i)) .eq. 0d0 ) then
        fwx(i) = 0d0
        fwy(i) = 0d0
    	mw(i) = - 2d0 * pi / 5d0 * r(i) ** 5 * rhowater * Cswater * &
				omega(i) * ABS(omega(i))
	else
		mw(i) = - 3d0 / 4d0 * pi * rhowater * Cswater * &
				sqrt( uw ** 2 + vw ** 2) * omega(i) * r(i) ** 4
	end if

end subroutine forcing


subroutine sheltering (j, i)

    implicit none

    include "parameter.h"
    include "CB_forcings.h"
    include "CB_variables.h"

    integer, intent(in) :: i, j

    double precision :: S_shelter

    ! sheltering height from air and water
    hsfa(j, i) = S_shelter(hfa(i), hfa(j), deltan(j,i), cosa(j,i), &
                    sina(j,i), ua, va)
    hsfw(j, i) = S_shelter(hfw(i), hfw(j), deltan(j,i), cosa(j,i), &
                    sina(j,i), uw, vw)

    ! sheltering for the reverse direction
    hsfa(j, i) = S_shelter(hfa(j), hfa(i), deltan(j,i), -cosa(j,i), &
                    -sina(j,i), ua, va)
    hsfw(j, i) = S_shelter(hfw(j), hfw(i), deltan(j,i), -cosa(j,i), &
                    -sina(j,i), uw, vw)
end subroutine


subroutine coriolis (i)

	implicit none

	integer, intent(in) :: i
    double precision :: omega_earth

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    omega_earth = 7.2921d-5

	fcorx(i)  = mass(i) * 2 * omega_earth * sin(lat) * v(i)

    fcory(i)  = - mass(i) * 2 * omega_earth * sin(lat) * u(i)

end subroutine coriolis


double precision function log_profile (hf, z0)

    implicit none

    double precision, intent(in) :: hf, z0

    log_profile = ( log( hf / z0 ) ** 2.0d0 - 2.0d0 *           &
                         log( hf / z0 ) + 2.0d0 - 2.0d0 * z0 / hf) / &
                         log( 10 / z0 ) ** 2.0d0

end function log_profile

double precision function heaviside (x)

    implicit none

    double precision, intent(in) :: x

    heaviside = sign(5d-1, x) + 0.5

end function heaviside

double precision function S_shelter (hfi, hfj, deltan, cosa, sina, uf, vf)

    implicit none

    double precision, intent(in) :: hfi, hfj, deltan
    double precision, intent(in) :: uf, vf, cosa, sina
    
    double precision :: heaviside

    S_shelter = (                                               &
                heaviside(deltan) * heaviside(hfi - hfj) *      &
                (1 - hfj / hfi) + ( 1 - heaviside(deltan) ) *   &
                ( 1 - exp(- 0.18 * abs(deltan) / hfj) )         &
                ) * (-uf * cosa - vf * sina) /                  &
                sqrt(uf**2 + vf**2) *                           &
                (                                               &
                heaviside(-cosa * uf) + heaviside(-sina * vf)   &
                )

end function S_shelter


double precision function L2norm (u, v)

    implicit none
    
    double precision, intent(in) :: u, v

    L2norm = sqrt( u ** 2d0 + v ** 2d0 )

end function L2norm
