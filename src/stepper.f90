subroutine stepper (tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j
    integer, intent(in) :: tstep
    double precision :: cosa, sina, ndot, sdot

    ! Before new time step following arrays which play transitional
    !role between global and incremental arrays shall be set to zero
    s   =  0d0
    sx  =  0d0
    sy  =  0d0
    w   =  0d0
    wx  =  0d0
    wy  =  0d0
    d   =  0d0
    dx  =  0d0
    dy  =  0d0
    l   =  0d0
    lx  =  0d0
    ly  =  0d0

    do i = 1, n
        do j = 1, n
            
            if ( j /= i ) then
                if (abs(                                 &
                        sqrt(                            &
                            ( xy(i,1) - xy(j,1) ) ** 2 + &
                            ( xy(i,2) - xy(j,2) ) ** 2   &
                            )                            &
                        ) < xy(i,3) + xy(j,3)            &
                    ) then

                    ! Components of unit vector ei=(cosa,sina) are:
                    cosa = ( xy(j,1) - xy(i,1) ) /            &
                           ( sqrt(                            &
                                 ( xy(i,1) - xy(j,1) ) ** 2 + &
                                 ( xy(i,2) - xy(j,2) ) ** 2   &
                                 )                            &
                           )
                    
                    sina = ( xy(j,2) - xy(i,2) ) /            &
                           ( sqrt(                            &
                                 ( xy(i,1) - xy(j,1) ) ** 2 + &
                                 ( xy(i,2) - xy(j,2) ) ** 2   &
                                 )                            &
                           )

                    ! Normal components of the relative velocities:
                    ndot = ( vu(i,1) - vu(j,1) ) * cosa + &
                           ( vu(i,2) - vu(j,2) ) * sina

                    ! Tangential components of the relative velocities:
                    sdot = ( ( vu(i,1) - vu(j,1) ) * sina +         &
                             ( vu(i,2) - vu(j,2) ) * ( - cosa ) ) - &
                           ( ome(i) * xy(i,3) + ome(j) * xy(j,3) )

                    ! Normal component of the damping force at contacts, d(j,i), force induced from
                    ! particle B to particle A at their contact point
                    d(j,i)  =  ndot * eta

                    dx(j,i) =  d(j,i) * ( - cosa )
                    dy(j,i) =  d(j,i) * ( - sina )

                    ! Shear component of damping force at contact point
                    l(j,i)  =  - 1 * ( sdot * eta2 )

                    ! Global components of L(j,i) in cartesian 
                    lx(j,i)  =  l(j,i) * sina
                    ly(j,i)  =  l(j,i) * ( - cosa )

                    ! Increments of the normal force
                    s(j,i)   =  ndot * dt * k

                    ! Global components of s(j,i) in cartesian
                    sx(j,i)  =  s(j,i) * ( - cosa )
                    sy(j,i)  =  s(j,i) * ( - sina )

                    ! Increments of the tangential forces
                    w(j,i)   =  - 1 * ( sdot * dt * ks )

                    !Global components of w(j,i) in cartesian
                    wx(j,i)  =  w(j,i) * sina
                    wy(j,i)  =  w(j,i) * ( - cosa )

                else
                    ! If 2 particles are separated from each other no contact force exist
                    s(j,i)      =  0d0
                    fsx(j,i)    =  0d0
                    fsy(j,i)    =  0d0
                    fw(j,i)     =  0d0
                    w(i,j)      =  0d0
                    wx(i,j)     =  0d0
                    wy(i,j)     =  0d0
                    fwx(j,i)    =  0d0
                    fwy(j,i)    =  0d0
                    d(j,i)      =  0d0
                    dx(j,i)     =  0d0
                    dy(j,i)     =  0d0
                    l(j,i)      =  0d0
                    lx(j,i)     =  0d0
                    ly(j,i)     =  0d0

                end if

            else
                s(j,i)      =  0d0
                fsx(j,i)    =  0d0
                fsy(j,i)    =  0d0
                fw(j,i)     =  0d0
                w(i,j)      =  0d0
                wx(i,j)     =  0d0
                wy(i,j)     =  0d0
                fwx(j,i)    =  0d0
                fwy(j,i)    =  0d0
                d(j,i)      =  0d0
                dx(j,i)     =  0d0
                dy(j,i)     =  0d0
                l(j,i)      =  0d0
                lx(j,i)     =  0d0
                ly(j,i)     =  0d0

            end if

        end do
    end do

    call contact_forces
    call viscous_terms

    call forcing (tstep)

    call euler

end subroutine stepper

subroutine euler

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! Update new position/angle of particles
    do i = 1, N

        xy(i,1) = xy(i,1) + vu(i,1) * dt
        xy(i,2) = xy(i,2) + vu(i,2) * dt
        teta(i) = teta(i) + ome(i) * dt

    end do

end subroutine euler