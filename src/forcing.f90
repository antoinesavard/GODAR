subroutine forcing (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    double precision :: uw, vw, ua, va
    double precision :: fdax(n), fday(n), fsax(i), fsay(i)
    double precision :: fdwx(n), fdwy(n), fswx(i), fswy(i)

    double precision :: log_wind_profile, H_shelter
    double precision :: hsfa, hsfw
    
    uw = 0d0
    vw = 0d0
    ua = 20d0
    va = 20d0

    ! sheltering height from air and water
    hsfa = H_shelter(hfa(i), -deltan(i,j), cosa(i,j), sina(i,j), 	&
						ua - u(i), va - v(i))
    hsfw = H_shelter(hfw(i), -deltan(i,j), cosa(i,j), sina(i,j), 	&
                        uw - u(i), vw - v(i))
    
    ! wind drag forcing
    fdax(i) = rhoair * Cdair * (hfa(i) - hsfa) * r(i) * (ua - u(i)) * &
                ABS(ua - u(i)) * log_wind_profile(hfa(i),             &
                max(z0w, hsfa))

    fday(i) = rhoair * Cdair * (hfa(i) - hsfa) * r(i) * (va - v(i)) * &
                ABS(va - v(i)) * log_wind_profile(hfa(i),             &
                max(z0w, hsfa))

    ! wind skin forcing
    fsax(i) = 0.5 * rhoair * pi * r(i) ** 2 * Csair * &
                (ua - u(i)) * ABS(ua - u(i))

    fsay(i) = 0.5 * rhoair * pi * r(i) ** 2 * Csair * &
                (va - v(i)) * ABS(va - v(i))

    ! water drag forcing
    fdwx(i) = rhowater * Cdwater * (hfw(i) - hsfw) * r(i) * (uw - &
                u(i)) * ABS(uw - u(i)) * log_wind_profile(hfw(i), &
                max(z0w, hsfw))

    fdwy(i) = rhowater * Cdwater * (hfw(i) - hsfw) * r(i) * (vw - &
                v(i)) * ABS(vw - v(i)) * log_wind_profile(hfw(i), &
                max(z0w, hsfw))

    ! water skin forcing
    fswx(i) = 0.5 * rhowater * pi * r(i) ** 2 * Cswater * &
                (uw - u(i)) * ABS(uw - u(i))

    fswy(i) = 0.5 * rhowater * pi * r(i) ** 2 * Cswater * &
                (vw - v(i)) * ABS(vw - v(i))

    ! total forcing from air and water
    fax(i) = fdax(i) + fsax(i)
    fay(i) = fday(i) + fsay(i)

    fwx(i) = fdwx(i) + fswx(i)
    fwy(i) = fdwy(i) + fswy(i)

    ! torque induced drag due to rotation of floes when no speed
	! if speed, use second expression valid for |U| >> |omega*r|
	if ( abs(ua) + abs(va) .eq. 0 ) then
		ma(i) = - pi / 5 * r(i) ** 5 * rhoair * Csair * omega(i) * &
				ABS(omega(i))
	else
		ma(i) = - 3 / 8 * pi * rhoair * Csair * &
				sqrt( ua ** 2 + va ** 2) * omega(i) * r ** 4
	end if

	if ( abs(uw) + abs(vw) .eq. 0 ) then
    	mw(i) = - pi / 5 * r(i) ** 5 * rhowater * Cswater * omega(i) * &
				ABS(omega(i))
	else
		ma(i) = - 3 / 8 * pi * rhowater * Cswater * &
				sqrt( uw ** 2 + vw ** 2) * omega(i) * r ** 4
	end if

end subroutine forcing


subroutine coriolis (i, j)

	implicit none

	integer, intent(in) :: i, j

end subroutine coriolis


double precision function log_wind_profile (hf, z0)

    implicit none

    double precision, intent(in) :: hf, z0

    log_wind_profile = ( log( hf / z0 ) ** 2.0d0 - 2.0d0 *           &
                         log( hf / z0 ) + 2.0d0 - 2.0d0 * z0 / hf) / &
                         log( 10 / z0 ) ** 2.0d0

end function log_wind_profile

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
                ) * ( 1 - sina ) * heaviside(sign(1d0, cosa)) 

end function H_shelter