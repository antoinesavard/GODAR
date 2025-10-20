subroutine verify_bc (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: i
    
    integer :: tracker1, tracker2

    tracker1 = 0
    tracker2 = 0

!-----------------------------------------
!   ooo
!   xxx
!-----------------------------------------         

    if ( y(i) - r(i) < 0 ) then

        tracker1 = 1

        call check_wall(i ,0 ,0 ,1)

    end if

!------------------------------------------
!   xxx
!   ooo
!------------------------------------------

    if ( y(i) + r(i) > ny ) then

        tracker1 = 1

        call check_wall(i ,0 ,1 ,1)

    end if

!------------------------------------------
!   xoo
!   xoo
!   xoo
!------------------------------------------

    ! but we also need the regular version for the cases where 
    ! particles are just left-right boundaries, but not in corners
    ! note that we do not need the nested ifs here, because it has
    ! been dealed with above.
    if ( x(i) - r(i) < 0 ) then

        tracker2 = 1

        call check_wall(i ,1 ,0 ,2)

    end if

!------------------------------------------
!   oox
!   oox
!   oox
!------------------------------------------

    if ( x(i) + r(i) > nx ) then

        tracker2 = 1

        call check_wall(i ,1 ,1 ,2)

    end if

!------------------------------------------
!   resets
!------------------------------------------

    if (tracker1 + tracker2 .eq. 0) then

        call reset_boundary (i, 0)

    end if

    if (tracker2 .eq. 0) then

        call reset_boundary (i, 2)

    end if

    if (tracker1 .eq. 0) then

        call reset_boundary (i, 1)

    end if

end subroutine verify_bc


subroutine check_wall(i, dir1, dir2, bd)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: i, dir1, dir2, bd

    double precision :: cosa_bc, sina_bc

    call contact_bc(i, dir1, dir2, bd)

    ! angles
    cosa_bc = dir1 * (dir2 - (1 - dir2))
    sina_bc = (1 - dir1) * (dir2 - (1 - dir2))

    ! update the forces applied by the boundaries on each particle
    fx_bc(i) = fx_bc(i) - fn_bc(i) * cosa_bc +    &
                          fr_bc(i) * sina_bc
    fy_bc(i) = fy_bc(i) - fn_bc(i) * sina_bc -    &
                          fr_bc(i) * cosa_bc

    ! update the moment applied by the boundaries on each particle
    m_bc(i) = m_bc(i) - mc_bc(i) - r(i) * ft_bc(i)

    ! compute the stress using cauchy stress formula due to the boundaries, always negative
    ! off diag are always 0
    sigxx_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * cosa_bc ** 2
    sigyy_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * sina_bc ** 2
    sigxy_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * cosa_bc * sina_bc
    sigyx_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * sina_bc * cosa_bc

    ! compute the pressure
    ! total contact area
    ta_bc(i) = ta_bc(i) + delt_ridge_bc(i) * h(i)
    
    ! pressure from contacts and bonds
    p_bc(i) = p_bc(i) - fn_bc(i) * delt_ridge_bc(i) * h(i)

end subroutine check_wall


