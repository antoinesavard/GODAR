subroutine bond_forces_surface (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

	integer, intent(in) :: i, j

    double precision :: gamnb, gamtb, gambb
    double precision :: mbending, mrolling
    double precision :: krb, gamrb
    double precision :: m_redu, i_redu, hmin
    double precision :: knb_eff, ktb_eff

    ! Relative displacements in bond frame
    deltanb(j,i) = veln(j,i) * dt + deltanb(j,i)
    deltatb(j,i) = velt(j,i) * dt + deltatb(j,i)

    ! relative angle for bending and twisting
    thetarelb(j,i) = omegarel(j,i) * dt + thetarelb(j,i)

    ! rolling stiffness due to bond
    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )
    i_redu =  inertia(i) * inertia(j) / ( inertia(i) + inertia(j) )
    hmin   =  min(h(i), h(j))

    ! effective stiffnesses for the bond forces
    knb_eff = (1d0 - damageb(j, i)) * knb(j, i)
    ktb_eff = (1d0 - damageb(j, i)) * ktb(j, i)

    ! effective viscosity
    gamnb = sqrt( 4d0 * knb_eff * sb(j,i) * m_redu )
    gamtb = sqrt( 4d0 * ktb_eff * sb(j,i) * m_redu )
    gambb = sqrt( 4d0 * ktb_eff * ib(j,i) * i_redu )

    ! forces are computed from linear elastic material law
    ! F = -kx-cu but x>0 is elongation so the force is supposed
    ! to bring back the particles towards equilibrium, so that
    ! the force must be positive too (F>0) for particle i (which)
    ! is the one on which we are centered. And the reverse for
    ! particle j (F<0). But we had a sign in stepper so that
    ! the signs are all gucci (F=kx+cu).
    fbn(j, i) = knb_eff * sb(j, i) * deltanb(j,i) &
                + gamnb * veln(j,i)
    fbt(j, i) = ktb_eff * sb(j, i) * deltatb(j,i) &
                + gamtb * velt(j,i)

	! moments for bending and twisting motion
    mbending = ktb_eff * ib(j, i) * thetarelb(j,i) &
                + gambb * omegarel(j,i)

    ! moments due to rolling
    mrolling = krb * thetarelb(j, i) + gamrb * omegarel(j, i)

    ! ensures no rolling if moment is too big
    if ( abs(thetarelb(j, i)) > (sqrt(3d0) * sigmacb_crit * &
         hb(j,i) + abs(fbn(j,i))) / knb_eff / rb(j,i) ) then
            
        mrolling = 0
        !mrolling = abs(fbn(j,i)) * rb(j,i) / 3 * &
                    !sign(1d0, omegarel(j,i))

    end if

    ! total moment due to bonds
    mbb(j, i) = mbending + mrolling

end subroutine bond_forces_surface


subroutine bond_forces_euler (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: i, j

    double precision :: EA, EI, L
    double precision :: eta_i, eta_a
    double precision :: k_axial, k_shear1, k_shear2, &
                        k_rot1, k_rot2
    double precision :: gam_axial, gam_shear1, gam_shear2, &
                        gam_rot1, gam_rot2
    double precision :: m_redu, i_redu

    ! Relative displacements in bond frame
    deltanb(j,i) = veln(j,i) * dt + deltanb(j,i)
    deltatb(j,i) = veltb(j,i) * dt + deltatb(j,i)

    ! angle relative to beam axis for bending and twisting
    thetarelb(j,i) = omega(i) * dt + thetarelb(j,i)
    thetarelb(i,j) = omega(j) * dt + thetarelb(i,j)

    ! reduced variables for the viscosity
    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )
    i_redu =  inertia(i) * inertia(j) / ( inertia(i) + inertia(j) )

    ! Beam length
    L = lb(j,i)

    ! Axial and bending rigidities with damage
    EA = (1d0 - damageb(j,i)) * eb * sb(j,i)
    EI = (1d0 - damageb(j,i)) * eb * ib(j,i)

    ! Euler–Bernoulli stiffness coefficients
    k_axial  = EA / L
    k_shear1 = 12d0 * EI / L**3
    k_shear2 = 6d0  * EI / L**2
    k_rot1   = 4d0  * EI / L
    k_rot2   = 2d0  * EI / L

    ! Euler–Bernoulli viscosity coefficients
    ! There is a choice to be made about where to put damage in the 
    ! viscosity, we choose to put it outside the square root so that 
    ! the relaxation time is preserved for all damage levels (E/\eta), 
    ! but it could be inside too.
    ! If the viscosity is outside the square root, then we preserve the 
    ! ratio between the stiffness and viscosity as the bond is damaged
    ! (the relaxation time is preserved for all damage levels), but if 
    ! the viscosity is inside the square root, then the viscosity 
    ! decreases faster than the stiffness as the bond is damaged, and 
    ! the relaxation time decreases with damage, which may be more 
    ! physical, but may lead to more instability in the numerical 
    ! scheme.
    eta_i = 2d0 * sqrt( EI * i_redu / L )
    eta_a = 2d0 * sqrt( EA * m_redu * L )
    gam_axial  = eta_a / L
    gam_shear1 = 12d0 * eta_i / L**3
    gam_shear2 = 6d0  * eta_i / L**2
    gam_rot1   = 4d0  * eta_i / L
    gam_rot2   = 2d0  * eta_i / L

    ! Axial force
    fbn(j,i) = k_axial * deltanb(j,i) + gam_axial * veln(j,i)

    ! Transverse shear force
    fbt(j,i) = k_shear1 * deltatb(j,i) &
                + k_shear2 * (thetarelb(j,i) + thetarelb(i,j)) &
                + gam_shear1 * veltb(j,i) &
                + gam_shear2 * (omega(i) + omega(j))

    ! Moments (Euler–Bernoulli)
    mbb(j,i) =  k_rot1 * thetarelb(j,i) + k_rot2 * thetarelb(i,j) &
                - k_shear2 * deltatb(j,i) &
                + gam_rot1 * omega(i) + gam_rot2 * omega(j) &
                - gam_shear2 * veltb(j,i)
                
    ! Newton's 3rd law for moments            
    mbb(i,j) =  k_rot2 * thetarelb(j,i) + k_rot1 * thetarelb(i,j) &
                - k_shear2 * deltatb(j,i) &
                + gam_rot2 * omega(i) + gam_rot1 * omega(j) &
                - gam_shear2 * veltb(j,i)

end subroutine bond_forces_euler


subroutine bond_breaking (j, i)

	implicit none

	include "parameter.h"
	include "CB_variables.h"
	include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: i, j
    double precision :: phi, psi

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

    psi = min( 1d0, 1d0 / sqrt( phi ) )

    ! compute damage in the bond
    damageb(j,i) = damageb(j,i) &
                    + (1d0 - psi) * (1d0 - damageb(j,i)) * dt / dtd &
                    - damageb(j,i) * dt / dth

    ! breaking of the bonds
    if ( damageb(j,i) .ge. dmax ) then

        bond(j, i) = 0
        fbn(j, i) = 0d0
        fbt(j, i) = 0d0
        mbb(j, i) = 0d0
        damageb(j, i) = 1d0

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
        damageb(j, i) = 0d0
        call bond_properties (j ,i)

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
    lb  (j, i) = dist(j, i)

    ! sb is the cross section area
    ! ib is the inertia
    sb  (j, i) = 2d0 * rb (j, i) * hb (j, i)
	ib  (j, i) = 2d0 / 3d0 * hb (j, i) * rb (j, i) ** 3d0

    ! knb is the normal k in Hooke's
    ! ktb is the tangent k in Hooke's
	! knb (j, i) = eb / lb (j, i)
	! ktb (j, i) = 5d0 / 6d0 * gb / lb (j, i)

end subroutine bond_properties