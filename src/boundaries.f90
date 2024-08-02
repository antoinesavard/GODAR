subroutine verify_bc (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    
    integer :: tracker1, tracker2
    double precision :: ft_bc_tmp, mc_bc_tmp

    tracker1 = 0
    tracker2 = 0

!-----------------------------------------
!   ooo
!   xxx
!-----------------------------------------         

    if ( y(i) - r(i) < 0 ) then

        tracker1 = 1

        call contact_bc (i, 0, 0, 1)

        ! update the forces applied by the boundaries on each particle
        fy_bc(i) = fn_bc(i)

        ! update the moment applied by the boundaries on each particle
        m_bc(i) = m_bc(i) + mc_bc(i) + r(i) * ft_bc(i)

    end if

!------------------------------------------
!   xxx
!   ooo
!------------------------------------------

    if ( y(i) + r(i) > ny ) then

        tracker1 = 1

        call contact_bc (i, 0, 1, 1)

        ! update the forces applied by the boundaries on each particle
        fy_bc(i) = -fn_bc(i)

        ! update the moment applied by the boundaries on each particle
        m_bc(i) = m_bc(i) + mc_bc(i) - r(i) * ft_bc(i)

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

        call contact_bc (i, 1, 0, 2)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = fn_bc(i)

        ! update the moment applied by the boundaries on each particle
        m_bc(i) = m_bc(i) + mc_bc(i) - r(i) * ft_bc(i) 

    end if

!------------------------------------------
!   oox
!   oox
!   oox
!------------------------------------------

    if ( x(i) + r(i) > nx ) then

        tracker2 = 1

        call contact_bc (i, 1, 1, 2)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = -fn_bc(i)

        ! update the moment applied by the boundaries on each particle
        m_bc(i) = m_bc(i) + mc_bc(i) + r(i) * ft_bc(i)

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