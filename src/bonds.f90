subroutine bond_forces (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

	integer, intent(in) :: i, j

    double precision :: mbending, mrolling
    double precision :: krb, knc
    double precision :: fit
    double precision :: r_redu, hmin, hb

    ! rolling stiffness due to bond
    r_redu =  r(i) * r(j) / ( r(i) + r(j) )

    hmin   =  min(h(i), h(j))

    knc    =  pi * ec * hmin *                  &
                fit( deltan(j,i) * r_redu /     &
                ( 2 * hmin ** 2 ) )

    krb    =  knc * deltat(j,i) ** 2 / 12

    ! forces are computed from linear elastic material law
    fbn(j, i) = -knb(j, i) * sb(j, i) * veln(j, i) * dt
    fbt(j, i) = -ktb(j, i) * sb(j, i) * velt(j, i) * dt

	! moments for bending and twisting motion
    mbending = -ktb(j, i) * ib(j, i) * omegarel(j, i) * dt

    ! moments due to rolling
    mrolling = -krb * omegarel(j, i) * dt

    ! ensures no rolling if moment is too big
    if ( abs(omegarel(j,i)) > 2 * (sqrt(3d0) * sigmacb_crit * hb(j,i) &
        + abs(fcn(j,i))) / knc * deltat(j,i) ) then
            
        mrolling = -abs(fcn(j,i)) * deltat(j,i) / 6 * &
                    sign(1d0, omegarel(j,i))

    end if

    ! total moment due to bonds
    mbb(j, i) = mbending + mrolling

end subroutine bond_forces


subroutine bond_breaking (j, i)

	implicit none

	include "parameter.h"
	include "CB_variables.h"
	include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: i, j

	taub(j, i) = abs(fbt(j, i)) / sb(j, i)
	sigmatb(j, i) = - fbn(j, i) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)
	sigmacb(j, i) = fbn(j, i) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)

	if ( taub(j, i) .gt. tau_crit * hb(j,i) ) then
		
		bond(j, i) = 0
        fbn(j, i)  = 0d0
        fbt(j, i)  = 0d0

	else if ( sigmacb(i, j) .gt. sigmacb_crit * hb(j,i) ) then
		
		bond(j, i) = 0
        fbn(j, i)  = 0d0
        fbt(j, i)  = 0d0

	else if ( sigmatb(j, i) .gt. sigmatb_crit * hb(j,i) ) then

		bond(j, i) = 0
        fbn(j, i)  = 0d0
        fbt(j, i)  = 0d0

	end if


end subroutine bond_breaking


subroutine bond_creation (j, i)

	implicit none

	include "parameter.h"
	include "CB_variables.h"
	include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: i, j

	if ( deltan(j, i) .ge. 0.01 * r(i)) then
		
		bond(j, i) = 1

	end if   

end subroutine bond_creation