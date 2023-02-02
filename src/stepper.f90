subroutine stepper (tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j
    integer, intent(in) :: tstep

    ! reinitialize force arrays
    do i = 1, n
        do j = 1, n
            tfx = 0d0
            tfy = 0d0
            fcn = 0d0
            fct = 0d0
        end do
    end do
    
    ! put yourself in the referential of the ith particle
    do i = 1, n - 1 
        do j = i + 1, n
    
            ! Components of unit vector ei=(cosa,sina) are:
            cosa(i,j) = ( x(j) - x(i) ) /                 	 &
                        ( sqrt(                              &
                                ( x(i) - x(j) ) ** 2 +       &
                                ( y(i) - y(j) ) ** 2         &
                                )                            &
                        )
            
            sina(i,j) = ( y(j) - y(i) ) /                 	 &
                        ( sqrt(                              &
                                ( x(i) - x(j) ) ** 2 +       &
                                ( y(i) - y(j) ) ** 2         &
                                )                            &
                        )

            ! Normal components of the relative velocities:
            veln(i,j) = ( u(j) - u(i) ) * cosa(i,j) +     &
                        ( v(j) - v(i) ) * sina(i,j)

            ! Tangential components of the relative velocities:
            velt(i,j) = ( u(j) - u(i) ) * sina(i,j) +   &
                        ( v(j) - v(i) ) * cosa(i,j) -   &
                        ( omega(i) * r(i) + omega(j) * r(j) )

            ! normal overlap (displacement) deltan >=0
            deltan(i,j)  =  r(i) + r(j) -          & 
                            sqrt(                  & 
                            ( x(i) - x(j) ) ** 2 + & 
                            ( y(i) - y(j) ) ** 2   & 
                            )

            ! verify if two particles are colliding
            if ( deltan(i,j) .ge. 0 ) then

                call contact_forces (i, j)

                tfx(i) = tfx(i) - fcn(i,j) * cosa(i,j)
                tfy(i) = tfy(i) - fcn(i,j) * sina(i,j)

            end if

            call forcing(i, j)

            call moment(i, j)

            tfx(i) = tfx(i) + fax(i) + fwx(i)
            tfy(i) = tfy(i) + fay(i) + fwy(i)

        end do
    end do

    call velocity

    call euler

end subroutine stepper

subroutine euler

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! Update new position/angle of particles
    do i = 1, n

        x(i) = x(i) + u(i) * dt
        y(i) = y(i) + v(i) * dt
        theta(i) = theta(i) + omega(i) * dt

        call bc_verify (i)

    end do

end subroutine euler
