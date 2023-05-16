subroutine contact_forces (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: i, j

    double precision :: m_redu, r_redu, hmin
    double precision :: fit
    double precision :: knc, ktc, gamn, gamt
    double precision :: krc, mrolling

    deltat(j,i) = sqrt( r(i) ** 2 - ( (dist(j,i) ** 2 - &
                    r(j) ** 2 + r(i) ** 2) / (2 * dist(j,i)) ) ** 2 )

    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )

    r_redu =  r(i) * r(j) / ( r(i) + r(j) )

    hmin   =  min(h(i), h(j))

    knc    = pi * ec * hmin  *                  &
                fit( deltan(j,i) * r_redu /     &
                ( 2 * hmin ** 2 ) )

    ktc    = 6d0 * gc / ec * knc

    krc    = knc * deltat(j,i) ** 2 / 12

    gamn   = -beta * sqrt( 5d0 * knc * m_redu )

    gamt   = -2d0 * beta * sqrt( 5d0 * gc / ec * knc * m_redu )

    ! compute the normal/tangent force
    ! normal force is F=kx-du it works because this is the force needed
    ! make an compression of x knowing that there is a damping that acts
    ! againts us, and F is the force that the particle j is feeling
    ! while -F is the force the particle i is feeling
    fcn(j,i) = knc * deltan(j,i) - gamn * veln(j,i)

    fct(j,i) = ktc * deltat(j,i) - gamt * velt(j,i)

	! make sure that disks are slipping if not enough normal force
    call coulomb (j, i)

    if ( bond (j, i) .eq. 0 ) then
        ! moments due to rolling
        mrolling = -krc * omegarel(j, i) * dt

        ! ensures no rolling if moment is too big
        if ( abs(omegarel(j,i)) > 2 * abs(fcn(j,i)) / knc * &
                deltat(j,i) ) then
                
            mrolling = -abs(fcn(j,i)) * deltat(j,i) / 6 * &
                        sign(1d0, omegarel(j,i))
        
        end if

        ! total moment due to rolling
        mcc(j, i) = mrolling 
    end if

end subroutine contact_forces


double precision function fit (xi)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    double precision, intent(in)  :: xi
    double precision :: p1, p2, p3, q1, q2

    p1 = 9117d-4
    p2 = 2722d-4
    p3 = 3324d-6
    q1 = 1524d-3
    q2 = 3159d-5

    fit = ( p1 * xi ** 2 + p2 * xi + p3 ) / ( xi ** 2 + q1 * xi + q2 )

end function fit
