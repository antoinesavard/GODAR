subroutine bond_forces (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

	integer, intent(in) :: i, j

    double precision :: mbending, mrolling
    double precision :: krb, gamrb
    double precision :: m_redu, r_redu, hmin

    thetarelb(j,i) = -omegarel(j,i) * dt + thetarelb(j,i)

    ! compression has delta>0
    deltanb(j,i) = -veln(j,i) * dt + deltanb(j,i)
    deltatb(j,i) = -velt(j,i) * dt + deltatb(j,i)

    ! rolling stiffness due to bond
    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )
    r_redu =  r(i) * r(j) / ( r(i) + r(j) )

    hmin   =  min(h(i), h(j))

    ! stiffness coefficient
    krb    =  knb(j,i) * rb(j,i) ** 2 / 3
    gamrb  = gamma_d * rb(j,i) ** 2 / 3

    ! forces are computed from linear elastic material law
    ! F = -kx-cu but x>0 is elongation so the force is supposed
    ! to bring back the particles towards equilibrium, so that
    ! the force must be positive too (F>0) for particle i (which)
    ! is the one on which we are centered. And the reverse for
    ! particle j (F<0). But we had a sign in stepper so that
    ! the signs are all gucci (F=kx+cu).
    fbn(j, i) = knb(j, i) * sb(j, i) * deltanb(j, i) &
                - gamma_d * veln(j,i)
    fbt(j, i) = ktb(j, i) * sb(j, i) * deltatb(j, i) &
                - gamma_d * velt(j,i)

	! moments for bending and twisting motion
    mbending = ktb(j, i) * ib(j, i) * thetarelb(j,i) &
                - gamma_d * omegarel(j,i)

    ! moments due to rolling
    mrolling = krb * thetarelb(j, i) - gamrb * omegarel(j, i)

    ! ensures no rolling if moment is too big
    if ( abs(thetarelb(j, i)) > (sqrt(3d0) * sigmacb_crit * &
         hb(j,i) + abs(fbn(j,i))) / knb(j,i) / rb(j,i) ) then
            
        mrolling = 0
        !mrolling = abs(fcn(j,i)) * deltat(j,i) / 6 * &
                    !sign(1d0, omegarel(j,i))

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
	sigmatb(j, i) = - abs(fbn(j, i)) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)
	sigmacb(j, i) = abs(fbn(j, i)) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)

	if ( taub(j, i) .gt. tau_crit * hb(j,i) ) then
		
        print*, "1"
		bond(j, i) = 0
        fbn(j, i)  = 0d0
        fbt(j, i)  = 0d0

	else if ( sigmacb(i, j) .gt. sigmacb_crit * hb(j,i) ) then
		
        print*, "2"
		bond(j, i) = 0
        fbn(j, i)  = 0d0
        fbt(j, i)  = 0d0

	else if ( sigmatb(j, i) .gt. sigmatb_crit * hb(j,i) ) then

        print*, "3"
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


subroutine bond_properties (j, i)

	implicit none

	include "parameter.h"
	include "CB_variables.h"
	include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: i, j
    double precision :: gb

    gb = eb / 2d0 / (1 + poiss_ratio)

    ! bond properties
    ! rb is "radius" such that 2rb is the width
    ! hb is the thickness
    ! lb is the lenght
	rb  (j, i) = lambda_rb * min(r(i), r(j))
	hb  (j, i) = (h(i) + h(j)) / 2d0
	lb  (j, i) = lambda_lb * (r(i) + r(j))

    ! sb is the cross section area
    ! ib is the inertia
    sb  (j, i) = 2d0 * rb (j, i) * hb (j, i)
	ib  (j, i) = 2d0 / 3d0 * hb (j, i) * rb (j, i) ** 3d0

    ! knb is the normal k in Hooke's
    ! ktb is the tangent k in Hooke's
	knb (j, i) = eb / lb (j, i)
	ktb (j, i) = 5d0 / 6d0 * gb / lb (j, i)

end subroutine bond_properties