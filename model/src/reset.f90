subroutine reset_contact (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: j, i

    thetarelc(j,i)  = 0d0
    deltat(j,i)     = 0d0
    delt_ridge(j,i) = 0d0

end subroutine reset_contact


subroutine reset_bond (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: j, i

    thetarelb(j,i) = 0d0
    deltanb(j,i)   = 0d0
    deltatb(j,i)   = 0d0

end subroutine reset_bond


subroutine reset_boundary (i, bd)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, bd

    if (bd .eq. 1) then

        theta_bc1(i)  = 0d0
        deltat_bc1(i) = 0d0
    
    else if (bd .eq. 2) then

        theta_bc2(i)  = 0d0
        deltat_bc2(i) = 0d0

    else    
    
        theta_bc1(i)  = 0d0
        theta_bc2(i)  = 0d0
        deltat_bc1(i) = 0d0
        deltat_bc2(i) = 0d0

    end if

    delt_ridge_bc(i) = 0d0
     

end subroutine reset_boundary



subroutine reset_forces

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"
    include "CB_diagnostics.h"

    integer :: i

    ! reinitialize force arrays for contact and bonds
    do i = 1, n
        mc(i)    = 0d0
        mb(i)    = 0d0
        fcx(i)   = 0d0
        fcy(i)   = 0d0
        fbx(i)   = 0d0
        fby(i)   = 0d0
        ! and total force arrays
        m(i)     = 0d0
		tfx(i)   = 0d0
        tfy(i)   = 0d0
        ! boundary forces
        mc_bc(i) = 0d0
        m_bc(i)  = 0d0
        fn_bc(i) = 0d0
        ft_bc(i) = 0d0
        fx_bc(i) = 0d0
        fy_bc(i) = 0d0
        ! stresses
        sigxx(i) = 0d0
        sigyy(i) = 0d0
        sigxy(i) = 0d0
        sigyx(i) = 0d0
        ! boundary stress
        sigxx_bc(i) = 0d0
        sigyy_bc(i) = 0d0
        sigxy_bc(i) = 0d0
        sigyx_bc(i) = 0d0
        ! pressure
        tp(i)  = 0d0
        tac(i) = 0d0
        tab(i) = 0d0
        pc(i)  = 0d0
        pb(i)  = 0d0
        ! boundary pressure
        p_bc(i)  = 0d0
        ta_bc(i) = 0d0

        ! mpi reduce
        ! contact
        fcx_r(i)  = 0d0
        fcy_r(i)  = 0d0
        mc_r(i)  = 0d0
        ! bond
        fbx_r(i)  = 0d0
        fby_r(i)  = 0d0
        mb_r(i)  = 0d0
        ! forcing
        fwx_r(i)  = 0d0
        fwy_r(i)  = 0d0
        mw_r(i)  = 0d0
        fax_r(i)  = 0d0
        fay_r(i)  = 0d0
        ma_r(i)  = 0d0
        fcorx_r(i)  = 0d0
        fcory_r(i)  = 0d0
        ! boundary
        fx_bc_r(i)  = 0d0
        fy_bc_r(i)  = 0d0
        m_bc_r(i)  = 0d0
        ! stress
        sigxx_r(i)  = 0d0
        sigyy_r(i)  = 0d0
        sigxy_r(i)  = 0d0
        sigyx_r(i)  = 0d0
        ! boundary stress
        sigxx_bc_r(i)  = 0d0
        sigyy_bc_r(i)  = 0d0
        sigxy_bc_r(i)  = 0d0
        sigyx_bc_r(i)  = 0d0
        ! pressure
        tac_r(i) = 0d0
        tab_r(i) = 0d0
        pc_r(i)  = 0d0
        pb_r(i)  = 0d0
        ! boundary pressure
        ta_bc_r(i) = 0d0
        p_bc_r(i)  = 0d0
        ! total arrays
        m_r(i)   = 0d0
		tfx_r(i) = 0d0
        tfy_r(i) = 0d0
        tsigxx_r(i) = 0d0
        tsigyy_r(i) = 0d0
        tsigxy_r(i) = 0d0
        tsigyx_r(i) = 0d0
        tp_r(i) = 0d0
    end do

end subroutine reset_forces


subroutine reset_shelter

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"

    ! reinitialize sheltering height
    hsfa = 1d0
    hsfw = 1d0
    
    ! mpi
    hsfa_r = 1d0
    hsfw_r = 1d0

end subroutine reset_shelter