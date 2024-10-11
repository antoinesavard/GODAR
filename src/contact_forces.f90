subroutine contact_forces (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_options.h"

    integer, intent(in) :: i, j

    double precision :: m_redu, r_redu, hmin
    double precision :: fit
    double precision :: knc, ktc, gamn, gamt
    double precision :: krc, mrolling

    deltat(j,i) = 2 * sqrt( r(i) ** 2 - ( (dist(j,i) ** 2 - &
                    r(j) ** 2 + r(i) ** 2) / (2 * dist(j,i)) ) ** 2 )

    thetarelc(j,i) = omegarel(j,i) * dt + thetarelc(j,i)

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
    fcn(j,i) = knc * deltan(j,i) + gamn * veln(j,i)

    fct(j,i) = ktc * deltat(j,i) + gamt * velt(j,i)

    ! verify if we are in the plastic case or not
    if ( ridging .eqv. .true. ) then
        if ( sigmanc_crit * hmin .le. fcn(j,i) / deltat(j,i) / hmin ) &
        then
            
            call plastic_contact (j, i, m_redu, hmin, krc)

        end if
    end if

	! make sure that disks are slipping if not enough normal force
    call coulomb (j, i)

    if ( bond (j, i) .eq. 0 ) then
        ! moments due to rolling
        mrolling = -krc * thetarelc(j, i) 

        ! ensures no rolling if moment is too big
        if ( abs( thetarelc(j, i) ) > 2 * abs(fcn(j,i)) / knc / &
            deltat(j,i) ) then
                
            mrolling = 0
            !mrolling = -abs(fcn(j,i)) * deltat(j,i) / 6 * &
                        !sign(1d0, omegarel(j,i))
        
        end if

        ! total moment due to rolling
        mcc(j, i) = mrolling 
    end if

end subroutine contact_forces


subroutine contact_bc (i, dir1, dir2, bd)

    ! This subroutine computes the contact forces between 
    ! the particles and the walls. It uses the same law
    ! as the one used for floe--floe interations. 
    ! 
    ! Arguments: 
    !   i    (int): particle id
    !   dir1 (int): vertical (0) or horizontal (1)
    !   dir2 (int): bottom-left (0) or top-right (1)
    !   bd   (int): one (1) or two (2) walls

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_options.h"

    integer, intent(in) :: i
    integer, intent(in) :: dir1, dir2, bd

    double precision :: fit
    double precision :: knc, ktc, gamn, gamt
    double precision :: krc
    double precision :: deltat_bc, deltan_bc
    double precision :: mrolling_bc

    ! compute the compression of the particle
    deltan_bc = ( (1 - dir2) * (r(i) - x(i)) +                     &
                        dir2 * (r(i) + x(i) - nx) ) * dir1 +       &
                ( (1 - dir2) * (r(i) - y(i)) +                     &
                        dir2 * (r(i) + y(i) - ny) ) * (1 - dir1)

    deltat_bc = sqrt( r(i) ** 2 - ( r(i) - deltan_bc ) ** 2 )

    ! update the relative angle that the particle makes with bd
    if (bd .eq. 1) then
        theta_bc1(i) = omega(i) * dt + theta_bc1(i)
    else if (bd .eq. 2) then
        theta_bc2(i) = omega(i) * dt + theta_bc2(i)
    end if

    ! compute the springs constant
    knc    = pi * ec * h(i)  *                  &
                fit( deltan_bc * r(i) /         &
                ( 2 * h(i) ** 2 ) )

    ktc    = 6d0 * gc / ec * knc

    krc    = knc * deltat_bc ** 2 / 12

    ! compute the dashpots constant
    gamn   = -beta * sqrt( 5d0 * knc * mass(i) )

    gamt   = -2d0 * beta * sqrt( 5d0 * gc / ec * knc * mass(i) )

    ! compute the normal/tangent force
    ! is using dir as a way to pick the proper velocity for 
    ! normal or tangent force.
    fn_bc(i) = knc * deltan_bc &
                + gamn * ( dir1 * u(i) + (1 - dir1) * v(i) )

    ! check if there is a parallel velocity, if not,
    ! no tangential force.
    if ( (1 - dir1) * u(i) + dir1 * v(i) .ne. 0 ) then
        ft_bc(i) = ktc * deltat_bc &
                    + gamt * ( (1 - dir1) * u(i) + dir1 * v(i) )
    end if

    ! verify if we are in the plastic case or not
    if ( ridging .eqv. .true. ) then
        if ( sigmanc_crit * h(i) .le. fn_bc(i) / deltat_bc / h(i) ) &
        then
            
            call plastic_contact_bc (i, deltan_bc, deltat_bc, krc, &
                                    dir1)

        end if
    end if

	! make sure that disks are slipping if not enough normal force
    call coulomb_bc (i, dir1)

    ! moments due to rolling
    if (bd .eq. 1) then

        mrolling_bc = -krc * theta_bc1(i) 

        ! ensures no rolling if moment is too big
        if ( abs( theta_bc1(i) ) > 2 * abs(fn_bc(i)) / knc / &
            deltat_bc ) then
                
            mrolling_bc = 0  
            !mrolling_bc = -abs(fn_bc(i)) * deltat_bc / 6 * &
                        !sign(1d0, omega(i))
        
        end if

    else if (bd .eq. 2) then

        mrolling_bc = -krc * theta_bc2(i)

        ! ensures no rolling if moment is too big
        if ( abs( theta_bc2(i) ) > 2 * abs(fn_bc(i)) / knc / &
            deltat_bc ) then
                
            mrolling_bc = 0 
            !mrolling_bc = -abs(fn_bc(i)) * deltat_bc / 6 * &
                        !sign(1d0, omega(i))
        
        end if

    else 

        print*, "Anormal boundary rolling scheme"
        mrolling_bc = 0d0

    end if

    ! total moment due to rolling
    mc_bc(i) = mrolling_bc 

end subroutine contact_bc


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