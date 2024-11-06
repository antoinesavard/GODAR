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

    ! compression has delta_t > 0
    deltat(j,i) = -velt(j,i) * dt + deltat(j,i)
    
    delt_ridge(j,i) = 2 * sqrt( r(i) ** 2 - ( (dist(j,i) ** 2 - &
                    r(j) ** 2 + r(i) ** 2) / (2 * dist(j,i)) ) ** 2 ) 

    ! relative angle
    thetarelc(j,i) = omegarel(j,i) * dt + thetarelc(j,i)

    ! reduced variables
    m_redu =  mass(i) * mass(j) / ( mass(i) + mass(j) )
    r_redu =  r(i) * r(j) / ( r(i) + r(j) )

    hmin   =  min(h(i), h(j))

    ! stiffness
    knc    = pi * ec * hmin  *                  &
                fit( deltan(j,i) * r_redu /     &
                ( 2 * hmin ** 2 ) )

    ktc    = 6d0 * gc / ec * knc

    krc    = knc * delt_ridge(j,i) ** 2 / 12

    ! damping
    gamn   = -beta * sqrt( 4d0 * knc * m_redu )

    gamt   = -beta * sqrt( 4d0 * ktc * m_redu )

    ! compute the normal/tangent force
    fcn(j,i) = knc * deltan(j,i) - gamn * veln(j,i)

    fct(j,i) = ktc * deltat(j,i) - gamt * velt(j,i)

    ! verify if we are in the plastic case or not
    if ( ridging .eqv. .true. ) then
        if ( sigmanc_crit * hmin .le. fcn(j,i) / delt_ridge(j,i) &
        / hmin ) then
            
            call plastic_contact (j, i, m_redu, hmin, ktc, krc, gamt)

        end if
    end if

	! make sure that disks are slipping if not enough normal force
    call coulomb (j, i, ktc, gamt)

    if ( bond (j, i) .eq. 0 ) then
        ! moments due to rolling
        mrolling = -krc * thetarelc(j, i) 

        ! ensures no rolling if moment is too big
        if ( abs( thetarelc(j, i) ) > 2 * abs(fcn(j,i)) / knc / &
            delt_ridge(j,i) ) then
                
            mrolling = 0d0
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
    !   dir2 (int): bottom/left (0) or top/right (1)
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
    double precision :: veln_bc, velt_bc
    double precision :: cosa_bc, sina_bc

    ! angles
    cosa_bc = dir1 * (dir2 - (1 - dir2))
    sina_bc = (1 - dir1) * (dir2 - (1 - dir2))

    ! Normal components of the relative velocities:
	veln_bc = -u(i) * cosa_bc - v(i) * sina_bc

    ! Tangential components of the relative velocities:
    velt_bc = u(i) * sina_bc - v(i) * cosa_bc - omega(i) * r(i)


    ! compute the compression of the particle
    deltan_bc = ( (1 - dir2) * (r(i) - x(i)) +                     &
                        dir2 * (r(i) + x(i) - nx) ) * dir1 +       &
                ( (1 - dir2) * (r(i) - y(i)) +                     &
                        dir2 * (r(i) + y(i) - ny) ) * (1 - dir1)

    deltat_bc = 0d0
    delt_ridge_bc(i) = sqrt( r(i) ** 2 - ( r(i) - deltan_bc ) ** 2 )

    ! compression has delta_t > 0
    if (bd .eq. 1) then
        deltat_bc1(i) = -velt_bc * dt + deltat_bc1(i)
        deltat_bc = deltat_bc1(i)
    else if (bd .eq. 2) then
        deltat_bc2(i) = -velt_bc * dt + deltat_bc2(i)
        deltat_bc = deltat_bc2(i)
    end if

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

    krc    = knc * delt_ridge_bc(i) ** 2 / 12

    ! compute the dashpots constant
    gamn   = -beta * sqrt( 4d0 * knc * mass(i) )

    gamt   = -beta * sqrt( 4d0 * ktc * mass(i) )

    ! compute the normal/tangent force
    ! is using dir as a way to pick the proper velocity for 
    ! normal or tangent force.
    fn_bc(i) = knc * deltan_bc - gamn * veln_bc

    ! tangential force
    ft_bc(i) = ktc * deltat_bc - gamt * velt_bc

    ! verify if we are in the plastic case or not
    if ( ridging .eqv. .true. ) then
        if ( sigmanc_crit * h(i) .le. fn_bc(i) / delt_ridge_bc(i) &
        / h(i) ) then
            
            call plastic_contact_bc (i, veln_bc, velt_bc, deltan_bc, &
                                        deltat_bc, delt_ridge_bc(i), &
                                        ktc, krc, gamt)

        end if
    end if

    ! moments due to rolling
    if (bd .eq. 1) then

        ! make sure that disks are slipping if not enough normal force
        call coulomb_bc (i, velt_bc, ktc, gamt, deltat_bc1(i))

        mrolling_bc = -krc * theta_bc1(i)

        ! ensures no rolling if moment is too big
        if ( abs( theta_bc1(i) ) > 2 * abs(fn_bc(i)) / knc / &
            delt_ridge_bc(i) ) then
                
            mrolling_bc = 0d0  
            !mrolling_bc = -abs(fn_bc(i)) * deltat_bc / 6 * &
                        !sign(1d0, omega(i))
        
        end if

    else if (bd .eq. 2) then

        ! make sure that disks are slipping if not enough normal force
        call coulomb_bc (i, velt_bc, ktc, gamt, deltat_bc2(i))

        mrolling_bc = -krc * theta_bc2(i)

        ! ensures no rolling if moment is too big
        if ( abs( theta_bc2(i) ) > 2 * abs(fn_bc(i)) / knc / &
            delt_ridge_bc(i) ) then
                
            mrolling_bc = 0d0
            !mrolling_bc = -abs(fn_bc(i)) * deltat_bc / 6 * &
                        !sign(1d0, omega(i))
        
        end if

    else 

        print*, "Abnormal boundary rolling scheme"
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