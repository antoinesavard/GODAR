subroutine broadcasting_ini (thread_requested, restart)
! this routine broadcasts the initialization variables
! to all the other mpi ranks

    use mpi_f08
    use mask_io, only: nx_mask, ny_mask, dx_mask, x_origin, y_origin, &
                       Lx, Ly, mask_proj, sdf

    use parameters
    use const
    use mpi_var
    use variables
    use bonds
    use options
    use forcings
    use diagnostics

    implicit none


    integer, intent(inout) :: thread_requested, restart
    
    !-------------------------------------------------------------------
    ! openmp variable
    !-------------------------------------------------------------------
    call mpi_bcast(thread_requested, 1, mpi_integer,         &
                    master, mpi_comm_world, ierr)
    
    ! set the restart variable
    call mpi_bcast(restart, 1, mpi_integer,             &
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
    call mpi_bcast(active, n, mpi_logical,              &
                    master, mpi_comm_world, ierr)

    !-------------------------------------------------------------------
    ! other variable broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(u, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(v, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(damageb, n * n, mpi_double_precision,&
                    master, mpi_comm_world, ierr)

    ! bond geometry and saved bending angles (in theta_offset until
    ! bond_restore), only read from the output files on restart
    if (restart .eq. 1) then
        call mpi_bcast(lb, n * n, mpi_double_precision,     &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(rb, n * n, mpi_double_precision,     &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(hb, n * n, mpi_double_precision,     &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(theta_offset, n * n, mpi_double_precision,   &
                        master, mpi_comm_world, ierr)
    end if

    !-------------------------------------------------------------------
    ! constants broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(mass, n, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(inertia, n, mpi_double_precision,    &
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
    call mpi_bcast(deltat, n * n, mpi_double_precision,           &
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
    ! individual wind and currents variables broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(ua_i, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(va_i, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(uw_i, n, mpi_double_precision,        &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(vw_i, n, mpi_double_precision,        &
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
    call mpi_bcast(flag_diag_stress, 1, mpi_logical,    &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(flag_diag_pressure, 1, mpi_logical,  &
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
    call mpi_bcast(sigmacb_crit, 1, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(tau_crit, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(bond_lim, 1, mpi_double_precision,       &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(dmax, 1, mpi_double_precision,           &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(dtd, 1, mpi_double_precision,            &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(dth, 1, mpi_double_precision,            &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(phi_int, 1, mpi_double_precision,        &
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

    !-------------------------------------------------------------------
    ! mask / signed-distance-field broadcast
    !-------------------------------------------------------------------
    call mpi_bcast(nx_mask, 1, mpi_integer,                 &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(ny_mask, 1, mpi_integer,                 &
                    master, mpi_comm_world, ierr)

    if (nx_mask > 0 .and. ny_mask > 0) then
        call mpi_bcast(dx_mask, 1, mpi_double_precision,    &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(x_origin, 1, mpi_double_precision,   &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(y_origin, 1, mpi_double_precision,   &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(Lx, 1, mpi_double_precision,         &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(Ly, 1, mpi_double_precision,         &
                        master, mpi_comm_world, ierr)
        call mpi_bcast(mask_proj, 8, mpi_character,         &
                        master, mpi_comm_world, ierr)

        if (.not. allocated(sdf)) then
            allocate(sdf(nx_mask, ny_mask))
        end if
        call mpi_bcast(sdf, nx_mask * ny_mask, mpi_real,    &
                        master, mpi_comm_world, ierr)
    end if

end subroutine broadcasting_ini


subroutine broadcast_shape

    ! this routine reduces the minimum sheltering heights and the
    ! overlap volume ridged during the step across all processes, so
    ! that every process applies the same shape changes to all the
    ! particles (apply_ridging) and the shapes stay identical everywhere

    use mpi_f08
    use mpi_counts_mod, only: counts

    use parameters
    use forcings
    use mpi_var
    use variables
    use options

    implicit none


    ! local variables
    double precision, allocatable :: recv_min_a(:), recv_min_w(:)

    if ( shelter .eqv. .true. ) then
        ! allocate recv buffers
        allocate(recv_min_a(local_n))
        allocate(recv_min_w(local_n))

        ! sendbuf: local_hsfa_min (n), recvbuf: recv_min_a (local_n)
        call mpi_reduce_scatter(local_hsfa_min, recv_min_a, counts, &
                mpi_double_precision, mpi_min, mpi_comm_world, ierr)

        call mpi_reduce_scatter(local_hsfw_min, recv_min_w, counts, &
                mpi_double_precision, mpi_min, mpi_comm_world, ierr)

        ! place received minima into hsfX_min_r at the local positions
        hsfa_min_r(local_disp+1 : local_disp + local_n) = recv_min_a
        hsfw_min_r(local_disp+1 : local_disp + local_n) = recv_min_w

        ! deallocate recv buffers
        deallocate(recv_min_a)
        deallocate(recv_min_w)
    end if

    ! overlap volume ridged by the pairs of all ranks
    if ( ridging .eqv. .true. ) then
        call mpi_allreduce(MPI_IN_PLACE, dvol, n,                   &
                mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    end if

end subroutine broadcast_shape


subroutine broadcast_total_forces

    ! this routine broadcasts (and reduce) the total forces and moments
    ! at the end of each time step so that each process can compute 
    ! their own time stepping using the total forces

    use mpi_f08

    use parameters
    use variables
    use const
    use bonds
    use forcings
    use options
    use mpi_var
    use diagnostics

    implicit none


    integer :: i, a, idx, field_num
    double precision, allocatable :: buf(:)

    ! 3 fields (total forces and moment), plus 4 stress and 1 pressure
    ! fields only when those diagnostics are on. Each rank filled its
    ! own particles only (the others are 0), so the sum gathers them
    field_num = 3
    if ( flag_diag_stress .eqv. .true. ) field_num = field_num + 4
    if ( flag_diag_pressure .eqv. .true. ) field_num = field_num + 1

    allocate(buf(field_num * n))

    ! all the fields of particle 1, then of particle 2, ...
    do i = 1, n
        idx = (i - 1) * field_num
        a = 1    ; buf(idx + a) = tfx_r(i)
        a = a + 1; buf(idx + a) = tfy_r(i)
        a = a + 1; buf(idx + a) = m_r(i)

        if ( flag_diag_stress .eqv. .true. ) then
            a = a + 1; buf(idx + a) = tsigxx_r(i)
            a = a + 1; buf(idx + a) = tsigyy_r(i)
            a = a + 1; buf(idx + a) = tsigxy_r(i)
            a = a + 1; buf(idx + a) = tsigyx_r(i)
        end if

        if ( flag_diag_pressure .eqv. .true. ) then
            a = a + 1; buf(idx + a) = tp_r(i)
        end if
    end do

    ! total forces, moments, stresses and pressure in one message
    call mpi_allreduce(MPI_IN_PLACE, buf, field_num * n,           &
            mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

    do i = 1, n
        idx = (i - 1) * field_num
        a = 1    ; tfx(i) = buf(idx + a)
        a = a + 1; tfy(i) = buf(idx + a)
        a = a + 1; m(i)   = buf(idx + a)

        if ( flag_diag_stress .eqv. .true. ) then
            a = a + 1; tsigxx(i) = buf(idx + a)
            a = a + 1; tsigyy(i) = buf(idx + a)
            a = a + 1; tsigxy(i) = buf(idx + a)
            a = a + 1; tsigyx(i) = buf(idx + a)
        end if

        if ( flag_diag_pressure .eqv. .true. ) then
            a = a + 1; tp(i) = buf(idx + a)
        end if
    end do

    deallocate(buf)

    ! particles that left the domain on any rank (remove_exited)
    call mpi_allreduce( &
    MPI_IN_PLACE, exited, n, mpi_logical, &
    mpi_lor, mpi_comm_world, ierr)


end subroutine broadcast_total_forces


subroutine force_reduction

    ! this routine reduces and broadcast to all processes the
    ! intermediate forces (contact, bond) and moments before each 
    ! process combines their section into total forecs and moments.
    ! this is needed because the program uses Newton's third law.

    use mpi_f08

    use parameters
    use variables
    use const
    use bonds
    use forcings
    use options
    use mpi_var
    use diagnostics

    implicit none


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

subroutine force_reduction_fast
    
    use mpi_f08
    use mpi_counts_mod

    use parameters
    use variables
    use const
    use bonds
    use forcings
    use options
    use mpi_var
    use diagnostics

    implicit none


    integer :: i, a, k, idx
    integer :: total_elems, my_recvcount, field_num
    integer, allocatable :: recvcounts(:)
    double precision, allocatable :: sendbuf(:), recvbuf(:)

    ! 17 force fields, plus 12 stress and 6 pressure fields only when
    ! those diagnostics are on
    field_num = 17
    if ( flag_diag_stress .eqv. .true. ) field_num = field_num + 12
    if ( flag_diag_pressure .eqv. .true. ) field_num = field_num + 6

    total_elems = field_num * n

    ! allocate buffers
    allocate(sendbuf(total_elems))
    allocate(recvcounts(n_ranks))

    do i = 1, n_ranks
        recvcounts(i) = counts(i) * field_num
    end do

    ! recvbuf size for this rank (only receive its local portion)
    my_recvcount = recvcounts(rank+1)
    allocate(recvbuf(my_recvcount))

    ! put all the fields one after the other
    !
    ! sendbuf = [fcx(1), fcy(1), ..., fcx(2), fcy(2), ..., fcx(n), ...]
    !            -------------------  -------------------  -----------
    !                particle 1           particle 2       particle n
    !
    do i = 1, n
        idx = (i - 1) * field_num
        a = 1    ; sendbuf(idx + a) = fcx(i)
        a = a + 1; sendbuf(idx + a) = fcy(i)
        a = a + 1; sendbuf(idx + a) = mc(i)

        a = a + 1; sendbuf(idx + a) = fbx(i)
        a = a + 1; sendbuf(idx + a) = fby(i)
        a = a + 1; sendbuf(idx + a) = mb(i)

        a = a + 1; sendbuf(idx + a) = fwx(i)
        a = a + 1; sendbuf(idx + a) = fwy(i)
        a = a + 1; sendbuf(idx + a) = mw(i)
        a = a + 1; sendbuf(idx + a) = fax(i)
        a = a + 1; sendbuf(idx + a) = fay(i)
        a = a + 1; sendbuf(idx + a) = ma(i)
        a = a + 1; sendbuf(idx + a) = fcorx(i)
        a = a + 1; sendbuf(idx + a) = fcory(i)

        a = a + 1; sendbuf(idx + a) = fx_bc(i)
        a = a + 1; sendbuf(idx + a) = fy_bc(i)
        a = a + 1; sendbuf(idx + a) = m_bc(i)

        if ( flag_diag_stress .eqv. .true. ) then
            a = a + 1; sendbuf(idx + a) = sigxx(i)
            a = a + 1; sendbuf(idx + a) = sigyy(i)
            a = a + 1; sendbuf(idx + a) = sigxy(i)
            a = a + 1; sendbuf(idx + a) = sigyx(i)

            a = a + 1; sendbuf(idx + a) = sigxx_bc(i)
            a = a + 1; sendbuf(idx + a) = sigyy_bc(i)
            a = a + 1; sendbuf(idx + a) = sigxy_bc(i)
            a = a + 1; sendbuf(idx + a) = sigyx_bc(i)

            a = a + 1; sendbuf(idx + a) = sigxx_aw(i)
            a = a + 1; sendbuf(idx + a) = sigyy_aw(i)
            a = a + 1; sendbuf(idx + a) = sigxy_aw(i)
            a = a + 1; sendbuf(idx + a) = sigyx_aw(i)
        end if

        if ( flag_diag_pressure .eqv. .true. ) then
            a = a + 1; sendbuf(idx + a) = tac(i)
            a = a + 1; sendbuf(idx + a) = tab(i)
            a = a + 1; sendbuf(idx + a) = pc(i)
            a = a + 1; sendbuf(idx + a) = pb(i)

            a = a + 1; sendbuf(idx + a) = ta_bc(i)
            a = a + 1; sendbuf(idx + a) = p_bc(i)
        end if
    end do

    ! each rank provides full sendbuf of size total_elems
    ! recvbuf gets only recvcounts(rank) elements
    call mpi_reduce_scatter( sendbuf, recvbuf, recvcounts, mpi_double_precision, mpi_sum, mpi_comm_world, ierr )

    ! unpack recvbuf into local *_r arrays for indices 
    ! first_iter:last_iter 
    ! recvbuf is in the same per-particle block order 
    ! only for the local indices in this rank.
    do i = 1, local_n
        idx = (i - 1) * field_num
        ! global particle index
        k = local_disp + i

        a = 1    ; fcx_r(k) = recvbuf(idx + a)
        a = a + 1; fcy_r(k) = recvbuf(idx + a)
        a = a + 1; mc_r(k)  = recvbuf(idx + a)

        a = a + 1; fbx_r(k) = recvbuf(idx + a)
        a = a + 1; fby_r(k) = recvbuf(idx + a)
        a = a + 1; mb_r(k)  = recvbuf(idx + a)

        a = a + 1; fwx_r(k) = recvbuf(idx + a)
        a = a + 1; fwy_r(k) = recvbuf(idx + a)
        a = a + 1; mw_r(k)  = recvbuf(idx + a)
        a = a + 1; fax_r(k) = recvbuf(idx + a)
        a = a + 1; fay_r(k) = recvbuf(idx + a)
        a = a + 1; ma_r(k)  = recvbuf(idx + a)
        a = a + 1; fcorx_r(k) = recvbuf(idx + a)
        a = a + 1; fcory_r(k) = recvbuf(idx + a)

        a = a + 1; fx_bc_r(k) = recvbuf(idx + a)
        a = a + 1; fy_bc_r(k) = recvbuf(idx + a)
        a = a + 1; m_bc_r(k)  = recvbuf(idx + a)

        if ( flag_diag_stress .eqv. .true. ) then
            a = a + 1; sigxx_r(k) = recvbuf(idx + a)
            a = a + 1; sigyy_r(k) = recvbuf(idx + a)
            a = a + 1; sigxy_r(k) = recvbuf(idx + a)
            a = a + 1; sigyx_r(k) = recvbuf(idx + a)

            a = a + 1; sigxx_bc_r(k) = recvbuf(idx + a)
            a = a + 1; sigyy_bc_r(k) = recvbuf(idx + a)
            a = a + 1; sigxy_bc_r(k) = recvbuf(idx + a)
            a = a + 1; sigyx_bc_r(k) = recvbuf(idx + a)

            a = a + 1; sigxx_aw_r(k) = recvbuf(idx + a)
            a = a + 1; sigyy_aw_r(k) = recvbuf(idx + a)
            a = a + 1; sigxy_aw_r(k) = recvbuf(idx + a)
            a = a + 1; sigyx_aw_r(k) = recvbuf(idx + a)
        end if

        if ( flag_diag_pressure .eqv. .true. ) then
            a = a + 1; tac_r(k) = recvbuf(idx + a)
            a = a + 1; tab_r(k) = recvbuf(idx + a)
            a = a + 1; pc_r(k)  = recvbuf(idx + a)
            a = a + 1; pb_r(k)  = recvbuf(idx + a)

            a = a + 1; ta_bc_r(k) = recvbuf(idx + a)
            a = a + 1; p_bc_r(k)  = recvbuf(idx + a)
        end if
    end do

    deallocate(sendbuf)
    deallocate(recvbuf)
    deallocate(recvcounts)

end subroutine force_reduction_fast


subroutine gather_bonds_to_master()

    use mpi_f08
    use mpi_counts_mod

    use parameters
    use mpi_var
    use bonds

    implicit none


    integer :: i, j, k, idx
    integer :: num_local_bonds, num_total_bonds
    integer, allocatable :: local_i(:), local_j(:)
    integer, allocatable :: bond_recvcounts(:), bond_displs(:)
    integer, allocatable :: all_i(:), all_j(:)
    integer, allocatable :: all_bond_counts(:)
    integer, allocatable :: state_recvcounts(:), state_displs(:)

    ! per-bond state, in this order: damageb, lb, rb, hb,
    ! alpha_total(j,i), theta_offset(j,i), theta_offset(i,j)
    integer, parameter :: nstate = 7
    double precision, allocatable :: local_state(:,:), all_state(:,:)

    !------------------------------------------------------------
    ! Count local bonds
    !------------------------------------------------------------
    num_local_bonds = 0
    do i = first_iter, last_iter
        do j = 1, n
            if (bond(j, i) == 1) then
                num_local_bonds = num_local_bonds + 1
            end if
        end do
    end do

    !------------------------------------------------------------
    ! Store local bond pairs
    !------------------------------------------------------------
    allocate(local_i(num_local_bonds), local_j(num_local_bonds))
    allocate(local_state(nstate, num_local_bonds))

    idx = 0
    do i = first_iter, last_iter
        do j = 1, n
            if (bond(j, i) == 1) then
                idx = idx + 1
                local_i(idx) = i
                local_j(idx) = j
                local_state(:, idx) = [ damageb(j, i), lb(j, i),      &
                                        rb(j, i), hb(j, i),           &
                                        alpha_total(j, i),            &
                                        theta_offset(j, i),           &
                                        theta_offset(i, j) ]
            end if
        end do
    end do

    !------------------------------------------------------------
    ! Gather counts on master
    !------------------------------------------------------------
    allocate(all_bond_counts(n_ranks))

    call mpi_gather(num_local_bonds, 1, mpi_integer,          &
                    all_bond_counts, 1, mpi_integer, 0,       &
                    mpi_comm_world, ierr)

    !------------------------------------------------------------
    ! Master prepares recvcounts & displs for Gatherv
    !------------------------------------------------------------
    if (rank .eq. master) then
        allocate(bond_recvcounts(n_ranks))
        allocate(bond_displs(n_ranks))

        bond_recvcounts = all_bond_counts
        bond_displs(1)  = 0

        do k = 2, n_ranks
            bond_displs(k) = bond_displs(k-1) + bond_recvcounts(k-1)
        end do

        num_total_bonds = bond_displs(n_ranks) + bond_recvcounts(n_ranks)

        allocate(all_i(num_total_bonds))
        allocate(all_j(num_total_bonds))
        allocate(all_state(nstate, num_total_bonds))

        ! the state buffers hold nstate values per bond
        allocate(state_recvcounts(n_ranks))
        allocate(state_displs(n_ranks))

        state_recvcounts = nstate * bond_recvcounts
        state_displs     = nstate * bond_displs
    end if

    !------------------------------------------------------------
    ! Gatherv the i-indices
    !------------------------------------------------------------
    call mpi_gatherv(local_i, num_local_bonds, mpi_integer,           &
                     all_i, bond_recvcounts, bond_displs, mpi_integer,&
                     0, mpi_comm_world, ierr)

    !------------------------------------------------------------
    ! Gatherv the j-indices
    !------------------------------------------------------------
    call mpi_gatherv(local_j, num_local_bonds, mpi_integer,           &
                     all_j, bond_recvcounts, bond_displs, mpi_integer,&
                     0, mpi_comm_world, ierr)

    !------------------------------------------------------------
    ! Gatherv the bond state values
    !------------------------------------------------------------
    call mpi_gatherv(local_state, nstate * num_local_bonds,           &
                     mpi_double_precision, all_state,                 &
                     state_recvcounts, state_displs,                  &
                     mpi_double_precision, 0, mpi_comm_world, ierr)

    !------------------------------------------------------------
    ! Master rank now has all_i(k), all_j(k) for k=1..num_total_bonds
    ! and can reconstruct the full bond and bond state arrays
    !------------------------------------------------------------
    if (rank .eq. master) then
        ! clear the columns owned by the other ranks, so that bonds
        ! broken there are removed before rebuilding from the list
        bond(:, 1:first_iter-1) = 0
        bond(:, last_iter+1:n)  = 0

        do k = 1, num_total_bonds
            i = all_i(k)
            j = all_j(k)

            bond(j, i)         = 1
            damageb(j, i)      = all_state(1, k)
            lb(j, i)           = all_state(2, k)
            rb(j, i)           = all_state(3, k)
            hb(j, i)           = all_state(4, k)
            alpha_total(j, i)  = all_state(5, k)
            theta_offset(j, i) = all_state(6, k)
            theta_offset(i, j) = all_state(7, k)
        end do
    end if

    ! deallocate
    deallocate(local_i, local_j)
    deallocate(local_state)
    if (rank .eq. master) then
        deallocate(bond_recvcounts, bond_displs)
        deallocate(state_recvcounts, state_displs)
        deallocate(all_i, all_j)
        deallocate(all_state)
    end if

end subroutine gather_bonds_to_master
