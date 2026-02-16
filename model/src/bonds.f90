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
    double precision :: knb_eff, ktb_eff

    thetarelb(j,i) = -omegarel(j,i) * dt + thetarelb(j,i)

    ! compression has delta>0
    deltanb(j,i) = -veln(j,i) * dt + deltanb(j,i)
    deltatb(j,i) = -velt(j,i) * dt + deltatb(j,i)

    ! rolling stiffness due to bond
    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )
    r_redu =  r(i) * r(j) / ( r(i) + r(j) )

    hmin   =  min(h(i), h(j))

    ! effective stiffnesses for the bond forces
    knb_eff = (1d0 - damageb(j, i)) * knb(j, i)
    ktb_eff = (1d0 - damageb(j, i)) * ktb(j, i)
    
    ! rolling stiffness coefficient
    krb    =  knb_eff * rb(j,i) ** 2 / 3
    gamrb  = gamma_d * rb(j,i) ** 2 / 3

    ! forces are computed from linear elastic material law
    ! F = -kx-cu but x>0 is elongation so the force is supposed
    ! to bring back the particles towards equilibrium, so that
    ! the force must be positive too (F>0) for particle i (which)
    ! is the one on which we are centered. And the reverse for
    ! particle j (F<0). But we had a sign in stepper so that
    ! the signs are all gucci (F=kx+cu).
    fbn(j, i) = knb_eff * sb(j, i) * deltanb(j, i) &
                - gamma_d * veln(j,i)
    fbt(j, i) = ktb_eff * sb(j, i) * deltatb(j, i) &
                - gamma_d * velt(j,i)

	! moments for bending and twisting motion
    mbending = ktb_eff * ib(j, i) * thetarelb(j,i) &
                - gamma_d * omegarel(j,i)

    ! moments due to rolling
    mrolling = krb * thetarelb(j, i) - gamrb * omegarel(j, i)

    ! ensures no rolling if moment is too big
    if ( abs(thetarelb(j, i)) > (sqrt(3d0) * sigmacb_crit * &
         hb(j,i) + abs(fbn(j,i))) / knb_eff / rb(j,i) ) then
            
        mrolling = 0
        !mrolling = abs(fbn(j,i)) * rb(j,i) / 3 * &
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
    double precision :: phi

    ! compute stresses in the bond
	taub(j, i) = abs(fbt(j, i)) / sb(j, i)
	sigmatb(j, i) = - fbn(j, i) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)
	sigmacb(j, i) = fbn(j, i) / sb(j, i) + abs(mbb(j, i)) * 	&
					rb(j, i) / ib(j, i)

    ! compute the failure criteria for the bond
	phi = (taub(j, i) / (tau_crit * hb(j,i))) ** 2d0

	if ( sigmacb(j, i) .gt. 0 ) then	
        phi = phi + (sigmacb(j, i) / (sigmacb_crit * hb(j,i))) ** 2d0
    end if

	if ( sigmatb(j, i) .gt. 0 ) then
        phi = phi + (sigmatb(j, i) / (sigmatb_crit * hb(j,i))) ** 2d0
	end if

    ! compute damage in the bond
    if (phi > 1d0) then
        damageb(j,i) = max(damageb(j,i), 1d0 - 1d0/phi)
    end if

    ! breaking of the bonds
    if ( damageb(j,i) .ge. 0.9d0 ) then

        bond(j, i) = 0
        fbn(j, i) = 0d0
        fbt(j, i) = 0d0
        damageb(j, i) = 1d0

    end if

    print*, "phi = ", phi, "sigmacb = ", sigmacb(j, i), "damage = ", damageb(j, i)


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
        damageb(j, i) = 0d0

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
	lb  (j, i) = lambda_lb * dist(j, i)

    ! sb is the cross section area
    ! ib is the inertia
    sb  (j, i) = 2d0 * rb (j, i) * hb (j, i)
	ib  (j, i) = 2d0 / 3d0 * hb (j, i) * rb (j, i) ** 3d0

    ! knb is the normal k in Hooke's
    ! ktb is the tangent k in Hooke's
	knb (j, i) = eb / lb (j, i)
	ktb (j, i) = 5d0 / 6d0 * gb / lb (j, i)

end subroutine bond_properties