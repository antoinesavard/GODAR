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
    ! been dealt with above.
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


subroutine verify_bc_mask (i)

    use mask_io, only: sdf_at, sdf_grad, sdf_outside

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: i

    double precision :: d, gx, gy, gnorm
    double precision :: cosa_bc, sina_bc, deltan_bc

    double precision, parameter :: gnorm_floor = 1d-3
    double precision, parameter :: outside_thr = 0.5d0 * sdf_outside

    if (.not. active(i)) return

    d = sdf_at(x(i), y(i))

    ! open-boundary exit: particle left the grid -> tag inactive,
    ! break all bonds involving i
    if (d > outside_thr) then
        active(i) = .false.
        bond(:, i) = 0
        bond(i, :) = 0
        call reset_boundary (i, 0)
        return
    end if

    if (d >= r(i)) then
        call reset_boundary (i, 0)
        return
    end if

    call sdf_grad(x(i), y(i), gx, gy)
    gnorm = sqrt(gx * gx + gy * gy)

    ! medial axis: gradient ill-defined; skip this step
    if (gnorm < gnorm_floor) then
        call reset_boundary (i, 0)
        return
    end if

    ! check_wall_mask / contact_bc_mask expect the normal pointing from
    ! the particle toward the wall; sdf gradient points the other way
    cosa_bc   = -gx / gnorm
    sina_bc   = -gy / gnorm
    deltan_bc = r(i) - d

    call check_wall_mask (i, cosa_bc, sina_bc, deltan_bc)

end subroutine verify_bc_mask


subroutine check_wall_mask (i, cosa_bc, sina_bc, deltan_bc)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: i
    double precision, intent(in) :: cosa_bc, sina_bc, deltan_bc

    call contact_bc_mask (i, cosa_bc, sina_bc, deltan_bc)

    fx_bc(i) = fx_bc(i) - fn_bc(i) * cosa_bc +    &
                          fr_bc(i) * sina_bc
    fy_bc(i) = fy_bc(i) - fn_bc(i) * sina_bc -    &
                          fr_bc(i) * cosa_bc

    m_bc(i) = m_bc(i) - mc_bc(i) - r(i) * ft_bc(i)

    sigxx_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * cosa_bc ** 2
    sigyy_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * sina_bc ** 2
    sigxy_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * cosa_bc * sina_bc
    sigyx_bc(i) = -sqrt(fn_bc(i) ** 2 + ft_bc(i) ** 2) * r(i) &
                    * sina_bc * cosa_bc

    ta_bc(i) = ta_bc(i) + delt_ridge_bc(i) * h(i)
    p_bc(i)  = p_bc(i)  - fn_bc(i) * delt_ridge_bc(i) * h(i)

end subroutine check_wall_mask


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


