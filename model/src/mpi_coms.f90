subroutine broadcasting_ini (num_threads)

    use mpi_f08

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_mpi.h"
    include "CB_variables.h"
    include "CB_bond.h"
    include "CB_options.h"
    include "CB_forcings.h"
    include "CB_diagnostics.h"

    integer, intent(inout) :: num_threads
    
    !-------------------------------------------------------------------
    ! openmp variable
    !-------------------------------------------------------------------
    call mpi_bcast(num_threads, 1, mpi_integer,         &
                    master, mpi_comm_world, ierr)

    !-------------------------------------------------------------------
    ! variable broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(x, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(y, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(r, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(h, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(theta, n, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(omega, n, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(bond, n * n, mpi_integer,            &
                    master, mpi_comm_world, ierr)

    !-------------------------------------------------------------------
    ! other variable broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(u, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(v, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)

    !-------------------------------------------------------------------
    ! constants broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(mass, n, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(hfa, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(hfw, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
                    
    !-------------------------------------------------------------------
    ! 2D variables broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(thetarelc, n * n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(thetarelb, n * n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(deltat, n * n, mpi_double_precision,           &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(deltanb, n * n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(deltatb, n * n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(hsfa, n * n, mpi_double_precision,             &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(hsfw, n * n, mpi_double_precision,             &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(ac, n * n, mpi_double_precision,               &
                    master, mpi_comm_world, ierr)      

    !-------------------------------------------------------------------
    ! boundary variables broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(theta_bc1, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(theta_bc2, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(deltat_bc1, n, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(deltat_bc2, n, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)

    !-------------------------------------------------------------------
    ! namelist variables broadcast
    !-------------------------------------------------------------------
    ! options
    call mpi_bcast(dynamics, 1, mpi_logical,            &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(slipping, 1, mpi_logical,            &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(thermodyn, 1, mpi_logical,           &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(cohesion, 1, mpi_logical,            &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(ridging, 1, mpi_logical,             &
                    master, mpi_comm_world, ierr) 
    call mpi_bcast(shelter, 1, mpi_logical,             &
                    master, mpi_comm_world, ierr)

    ! numerical_param
    call mpi_bcast(rtree, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(ntree, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(dt, 1, mpi_double_precision,         &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(nt, 1, mpi_double_precision,         &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(comp, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)   

    ! physical_param
    call mpi_bcast(Cdair, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(Csair, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(Cdwater, 1, mpi_double_precision,    &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(Cswater, 1, mpi_double_precision,    &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(z0w, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(lat, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(rhoair, 1, mpi_double_precision,     &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(rhoice, 1, mpi_double_precision,     &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(rhowater, 1, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)       

    ! disk_param
    call mpi_bcast(e_modul, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(poiss_ratio, 1, mpi_double_precision,    &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(friction_coeff, 1, mpi_double_precision, &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(rest_coeff, 1, mpi_double_precision,     &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(sigmanc_crit, 1, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)

    ! bond_param
    call mpi_bcast(eb, 1, mpi_double_precision,             &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(lambda_rb, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(lambda_lb, 1, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(sigmatb_crit, 1, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(sigmacb_crit, 1, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(tau_crit, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(gamma_d, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(bond_lim, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)

    ! forcings
    call mpi_bcast(uw, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(vw, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(ua, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(va, 1, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)

    call mpi_bcast(pfn, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(pfs, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)

    ! other variables computed from namelist values
    call mpi_bcast(t, 1, mpi_double_precision,              &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(ec, 1, mpi_double_precision,             &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(gc, 1, mpi_double_precision,             &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(beta, 1, mpi_double_precision,           &
                    master, mpi_comm_world, ierr)   

end subroutine broadcasting_ini


subroutine broadcast_shape

    use mpi_f08
    use mpi_counts_mod, only: counts, displs

    implicit none

    include "parameter.h"
    include "CB_forcings.h"
    include "CB_mpi.h"
    include "CB_variables.h"

    ! local variables for allgatherv
    double precision, allocatable :: local_h(:), local_r(:),   &
                                    local_hfa(:), local_hfw(:)

    ! allocations
    allocate(local_h(local_n))
    local_h = h(first_iter:last_iter)

    allocate(local_r(local_n))
    local_r = r(first_iter:last_iter)

    allocate(local_hfa(local_n))
    local_hfa = hfa(first_iter:last_iter)

    allocate(local_hfw(local_n))
    local_hfw = hfw(first_iter:last_iter)

    ! sheltering coefficient
    call mpi_allreduce( &
    hsfa, hsfa_r, n * n, mpi_double_precision, &
    mpi_prod, mpi_comm_world, ierr)

    call mpi_allreduce( &
    hsfw, hsfw_r, n * n, mpi_double_precision, &
    mpi_prod, mpi_comm_world, ierr)

    ! thickness and radius
    call mpi_allgatherv(h(first_iter:last_iter),          &
                        local_n, mpi_double_precision,    &
                        h, counts, displs, mpi_double_precision,  &
                        mpi_comm_world, ierr )

    call mpi_allgatherv(r(first_iter:last_iter),          &
                        local_n, mpi_double_precision,    &
                        r, counts, displs, mpi_double_precision,  &
                        mpi_comm_world, ierr )

    ! freeboard heights
    call mpi_allgatherv(hfa(first_iter:last_iter),        &
                        local_n, mpi_double_precision,    &
                        hfa, counts, displs, mpi_double_precision,  &
                        mpi_comm_world, ierr )

    call mpi_allgatherv(hfw(first_iter:last_iter),        &
                        local_n, mpi_double_precision,    &
                        hfw, counts, displs, mpi_double_precision,  &
                        mpi_comm_world, ierr )

    ! deallocation
    deallocate(local_h)
    deallocate(local_r)
    deallocate(local_hfa)
    deallocate(local_hfw)

end subroutine broadcast_shape


subroutine broadcast_total_forces

    ! this routine broadcasts (and reduce) the total forces and moments
    ! at the end of each time step so that each process can compute 
    ! their own time stepping using the total forces

    use mpi_f08

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"
    include "CB_diagnostics.h"

    ! total forces
    call mpi_allreduce( &
    tfx_r, tfx, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    tfy_r, tfy, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    m_r, m, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! total stresses
    call mpi_allreduce( &
    tsigxx_r, tsigxx, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce( &
    tsigyy_r, tsigyy, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    tsigxy_r, tsigxy, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    tsigyx_r, tsigyx, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! total pressure
    call mpi_allreduce( &
    tp_r, tp, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    

end subroutine broadcast_total_forces


subroutine force_reduction

    ! this routine reduces and broadcast to all processes the
    ! intermediate forces (contact, bond) and moments before each 
    ! process combines their section into total forecs and moments.
    ! this is needed because the program uses Newton's third law.

    use mpi_f08

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"
    include "CB_diagnostics.h"

    ! contact
    call mpi_allreduce( &
    fcx, fcx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fcy, fcy_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    mc, mc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! bond
    call mpi_allreduce( &
    fbx, fbx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fby, fby_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    mb, mb_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! forcing
    call mpi_allreduce( &
    fwx, fwx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fwy, fwy_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    mw, mw_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fax, fax_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fay, fay_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    ma, ma_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fcorx, fcorx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fcory, fcory_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! boundary
    call mpi_allreduce( &
    fx_bc, fx_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    fy_bc, fy_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    m_bc, m_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! stress
    call mpi_allreduce( &
    sigxx, sigxx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce( &
    sigyy, sigyy_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigxy, sigxy_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigyx, sigyx_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! boundary stress
    call mpi_allreduce( &
    sigxx_bc, sigxx_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce( &
    sigyy_bc, sigyy_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigxy_bc, sigxy_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigyx_bc, sigyx_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! forcing stress
    call mpi_allreduce( &
    sigxx_aw, sigxx_aw_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce( &
    sigyy_aw, sigyy_aw_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigxy_aw, sigxy_aw_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    sigyx_aw, sigyx_aw_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! pressure
    call mpi_allreduce( &
    tac, tac_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce( &
    tab, tab_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    pc, pc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    pb, pb_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    ! boundary pressure
    call mpi_allreduce( &
    ta_bc, ta_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)

    call mpi_allreduce( &
    p_bc, p_bc_r, n, mpi_double_precision, &
    mpi_sum, mpi_comm_world, ierr)


end subroutine force_reduction

