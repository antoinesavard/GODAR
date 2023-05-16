subroutine forcing (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_forcings.h"

    integer, intent(in) :: i

    double precision :: fdax, fday, fsax, fsay
    double precision :: fdwx, fdwy, fswx, fswy

    double precision :: log_profile
    double precision :: hsfa_f(n), hsfw_f(n)

    hsfa_f = maxval(hsfa, dim=1)
    hsfw_f = maxval(hsfw, dim=1)

    if ( hsfa_f(i) .gt. 1d0 ) then
        hsfa_f(i) = 1d0
    end if

    if ( hsfw_f(i) .gt. 1d0 ) then
        hsfw_f(i) = 1d0
    end if    
    
    ! wind drag forcing
    fdax     = rhoair * Cdair * hfa(i) * r(i) * (ua - u(i)) *    &
                ABS(ua - u(i)) * log_profile(hfa(i),             &
                max(z0w, hsfa_f(i))) * hsfa_f(i)

    fday     = rhoair * Cdair * hfa(i) * r(i) * (va - v(i)) *    &
                ABS(va - v(i)) * log_profile(hfa(i),             &
                max(z0w, hsfa_f(i))) * hsfa_f(i)

    ! wind skin forcing
    fsax     = 0.5 * rhoair * pi * r(i) ** 2 * Csair * &
                (ua - u(i)) * ABS(ua - u(i))

    fsay     = 0.5 * rhoair * pi * r(i) ** 2 * Csair * &
                (va - v(i)) * ABS(va - v(i))

    ! water drag forcing
    fdwx     = rhowater * Cdwater * hfw(i) * r(i) * (uw -            &
                u(i)) * ABS(uw - u(i)) * log_profile(hfw(i),         &
                max(z0w, hsfw_f(i))) * hsfw_f(i)

    fdwy     = rhowater * Cdwater * hfw(i) * r(i) * (vw -            &
                v(i)) * ABS(vw - v(i)) * log_profile(hfw(i),         &
                max(z0w, hsfw_f(i))) * hsfw_f(i)

    ! water skin forcing
    fswx     = 0.5 * rhowater * pi * r(i) ** 2 * Cswater * &
                (uw - u(i)) * ABS(uw - u(i))

    fswy     = 0.5 * rhowater * pi * r(i) ** 2 * Cswater * &
                (vw - v(i)) * ABS(vw - v(i))

    ! total forcing from air and water
    fax(i) = fdax + fsax
    fay(i) = fday + fsay

    fwx(i) = fdwx + fswx
    fwy(i) = fdwy + fswy

    ! torque induced drag due to rotation of floes when no speed
	! if speed, use second expression valid for |U| >> |omega*r|
	if ( abs(ua) + abs(va) .eq. 0 ) then
		ma(i) = - pi / 5d0 * r(i) ** 5 * rhoair * Csair * omega(i) * &
				ABS(omega(i))
	else
		ma(i) = - 3d0 / 8d0 * pi * rhoair * Csair * &
				sqrt( ua ** 2 + va ** 2) * omega(i) * r(i) ** 4
	end if

	if ( abs(uw) + abs(vw) .eq. 0 ) then
    	mw(i) = - pi / 5d0 * r(i) ** 5 * rhowater * Cswater * &
				omega(i) * ABS(omega(i))
	else
		mw(i) = - 3d0 / 8d0 * pi * rhowater * Cswater * &
				sqrt( uw ** 2 + vw ** 2) * omega(i) * r(i) ** 4
	end if

end subroutine forcing


subroutine sheltering (j, i)

    implicit none

    include "parameter.h"
    include "CB_forcings.h"
    include "CB_variables.h"

    integer, intent(in) :: i, j

    double precision :: H_shelter

    ! sheltering height from air and water in units of h(i)
    hsfa(j, i) = H_shelter(hfa(j), -deltan(j,i), cosa(j,i), sina(j,i), &
                    ua - u(i), va - v(i)) / hfa(i)
    hsfw(j, i) = H_shelter(hfw(j), -deltan(j,i), cosa(j,i), sina(j,i), &
                    uw - u(i), vw - v(i)) / hfw(i)

end subroutine


subroutine coriolis (i)

	implicit none

	integer, intent(in) :: i

    include "parameter.h"
    include "CB_variables.h"

	!fcorx(i)  = mass(i) * f * v(i)

    !fcory(i)  = - mass(i) * f * u(i)

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

double precision function H_shelter (hf, deltan, cosa, sina, uf, vf)

    implicit none

    double precision, intent(in) :: hf, deltan
    double precision, intent(in) :: uf, vf, cosa, sina
    
    double precision :: u
    double precision :: heaviside

    u = sqrt(uf ** 2 + vf ** 2)

    H_shelter = (                                               &
                hf * ( 1 - heaviside(deltan) ) +                &
                hf * heaviside(deltan) * ( 1 - exp(-1d0) ) **   &
                ( 5d1 / (u + 1d-5) * heaviside(25d0 - u) +      &
                2d0 * heaviside(u - 25d0) )                     &
                ) *                                             &
                (                                               &
                heaviside(sign(1d0, -cosa * sign(1d0, uf)))      &
                * ( 1 - abs(sina) ) + ( 1 - abs(cosa) ) *       &
                heaviside(sign(1d0, -sina * sign(1d0, vf)))      &
                )

end function H_shelter