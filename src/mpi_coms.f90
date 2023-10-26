subroutine broadcasting_ini (num_threads)

    use mpi

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_mpi.h"
    include "CB_variables.h"
    include "CB_bond.h"
    include "CB_options.h"


    integer, intent(in) :: num_threads
    
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
    ! namelist variables broadcast
    !-------------------------------------------------------------------
    ! options
    call mpi_bcast(dynamics, 1, mpi_logical,            &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(thermodyn, 1, mpi_logical,           &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(cohesion, 1, mpi_logical,            &
                    master, mpi_comm_world, ierr)   
    call mpi_bcast(ridging, 1, mpi_logical,             &
                    master, mpi_comm_world, ierr) 

    ! numerical_param
    call mpi_bcast(rtree, 1, mpi_double_precision,      &
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

end subroutine broadcasting_ini


subroutine gathering_tstep

    use mpi

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"

    call mpi_gatherv( &
    fcx(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcx, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fcy(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcy, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fbx(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fbx, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fby(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fby, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fax(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcx, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fay(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcy, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fwx(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcx, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fwy(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcy, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fcorx(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcx, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    fcory(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    fcy, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    mc(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    mc, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

    call mpi_gatherv( &
    mb(first_iter:last_iter), iter_per_rank, mpi_double_precision, &
    mb, counts, disp, mpi_double_precision, &
    master, mpi_comm_world, ierr)

end subroutine gathering_tstep


subroutine broadcast_forces

    use mpi

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"

    call mpi_bcast(tfx, n, mpi_double_precision,    &
                    master, mpi_comm_world, ierr)

    call mpi_bcast(tfy, n, mpi_double_precision,    &
                    master, mpi_comm_world, ierr)

    call mpi_bcast(m, n, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)

end subroutine broadcast_forces

