subroutine stepper (tstep, restart)

    use mpi_f08
    use omp_lib
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_strings, only: str
    use kdtree_utils, only: neighbours, build_neighbours, neighbours_outdated
    use mask_io, only: nx_mask

    use parameters
    use variables
    use const
    use bonds
    use forcings
    use options
    use mpi_var
    use diagnostics
    use pairs, only: pair_t
    use timers, only: timer_on, timer_off, t_integ, &
    t_tree,t_pairs, t_forcing, t_comm

    implicit none

    integer :: i, j, k
    integer, intent(in) :: tstep, restart

    ! values computed for the current pair (j,i)
    type(pair_t) :: p

    ! sheltering coeff per thread arrays, (particle, thread) so that
    ! each thread writes in its own column
    integer :: thread_num, thread_id
    double precision, allocatable :: local_hsfa_min_thread(:,:)
    double precision, allocatable :: local_hsfw_min_thread(:,:)

    ! allocate local sheltering reduction arrays
    thread_num = omp_get_max_threads()

    if ( shelter .eqv. .true. ) then
        allocate(local_hsfa_min_thread(n, thread_num))
        allocate(local_hsfw_min_thread(n, thread_num))

        local_hsfa_min_thread = 1.0d0
        local_hsfw_min_thread = 1.0d0
    end if

    ! Velocity Verlet: advance positions before force computation
    ! (tree and contacts will be evaluated at x^{n+1})
    if ( tstep .ge. 1 ) then
        call timer_on(t_integ)
        call position
        call timer_off(t_integ)
    end if

    ! neighbour lists: the particles j > i within r(i) + rtree, rebuilt
    ! every ntree steps, or earlier if a particle moved more than half of
    ! the margin rtree - max(r) (a pair not listed could then touch)
    call timer_on(t_tree)
    if ( mod(tstep, int(ntree)) == 0 .or. &
         neighbours_outdated(x, y, rtree - maxval(r)) ) then
        call build_neighbours(x, y, r + rtree, first_iter, last_iter)
    end if
    call timer_off(t_tree)

    call timer_on(t_pairs)

    ! reset the forces and sheltering height
    call reset_forces

    ! reset the sheltering height to the maximum value (i.e. no
    ! sheltering) before computing the minimum value from the neighbors
    call reset_shelter
    
    ! put yourself in the referential of the ith particle
	! loop through all j particles and compute interactions

    !$omp parallel &
    !$omp private(i,j, thread_id, p) &
    !$omp reduction(+:fcx,fcy,mc,fbx,fby,mb) &
    !&&#ifdef DIAG
    !$omp reduction(+:sigxx,sigyy,sigxy,sigyx) &
    !$omp reduction(+:tac,tab,pc,pb)
    !&&#endif
    thread_id = omp_get_thread_num() + 1
    !$omp do schedule(dynamic, 1)
    do i = first_iter, last_iter

        if (.not. active(i)) cycle

        ! calculate the winds and currents applied on particle i
        call winds_currents(i)

        ! loop over the neighbours j > i of particle i (lower triangular
        ! matrices only)
        do k = 1, size(neighbours(i)%j)
            j = neighbours(i)%j(k)

            if (.not. active(j)) cycle

			! compute relative position and velocity
            p = pair_t()
            call rel_pos_vel (j, i, p)

			! bond initialization
            if ( tstep .eq. 0 .and. restart .ne. 1 ) then
                if ( cohesion .eqv. .true. ) then
                    call bond_creation (j, i, p)
                end if
			end if

            ! bond restoration from the restart files
            if ( tstep .eq. 0 .and. restart .eq. 1 ) then
                if ( bond(j, i) .eq. 1 ) then
                    call bond_restore (j, i, p)
                end if
            end if

            ! verify if two particles are colliding
            if ( p%deltan .gt. 0 ) then

!               call dilation (j, i) ! to implement
             
                call contact_forces (j, i, p)
!               call bond_creation (j, i) ! to implement
                
                ! change coordinate system
				! update contact force on particle i by particle j
                fcx(i) = fcx(i) - p%fcn * p%cosa +    &
                                        p%fcr * p%sina
                fcy(i) = fcy(i) - p%fcn * p%sina -    &
                                        p%fcr * p%cosa

                ! update moment on particule i by particule j due to tangent contact
                mc(i) = mc(i) - r(i) * p%fct - p%mcc

                ! Newton's third law
                ! update contact force on particle j by particle i
                fcx(j) = fcx(j) + p%fcn * p%cosa -    &
                                        p%fcr * p%sina
                fcy(j) = fcy(j) + p%fcn * p%sina +    &
                                        p%fcr * p%cosa

                ! update moment on particule j by particule i due to tangent contact
                mc(j) = mc(j) - r(j) * p%fct + p%mcc

                if ( flag_diag_pressure .eqv. .true. ) then
                    ! compute the average pressure inside particle i
                    !---------------------------------------------------
                    !
                    ! P_i = \sum_{c} Fcn_{ij} * a_{ij} / \sum_{c} a_{ij}
                    !
                    !---------------------------------------------------
                    block
                        double precision :: ac_ij

                        ! local area
                        ac_ij = p%delt_ridge * min(h(i), h(j))
                        ! total contact area
                        tac(i)  = tac(i) + ac_ij
                        ! pressure from contacts
                        pc(i)   = pc(i) - p%fcn * ac_ij
                        ! symmetric part
                        tac(j) = tac(j) + ac_ij
                        pc(j)  = pc(j) - p%fcn * ac_ij
                    end block
                end if

            else
            
                call reset_contact (j, i)

            end if

			! compute forces from bonds between particle i and j
			if ( bond (j, i) .eq. 1 ) then

				call bond_forces_timoshenko (j, i, p)
				call bond_breaking (j, i, p)

                if ( bond (j, i) .eq. 1 ) then

                    ! change coordinate system
                    ! update force on particle i by j due to bond
                    fbx(i) = fbx(i) - p%fbn * p%cosa +    &
                                        p%fbt * p%sina
                    fby(i) = fby(i) - p%fbn * p%sina -    &
                                        p%fbt * p%cosa

                    ! update moment on particule i by j to to bond
                    mb(i) = mb(i) + p%mbb_ji

                    ! Newton's third law
                    ! update force on particle j by i due to bond
                    fbx(j) = fbx(j) + p%fbn * p%cosa -    &
                                        p%fbt * p%sina
                    fby(j) = fby(j) + p%fbn * p%sina +    &
                                        p%fbt * p%cosa


                    ! update moment on particule j by i due to bond
                    mb(j) = mb(j) + p%mbb_ij

                    if ( flag_diag_pressure .eqv. .true. ) then
                        ! compute the average pressure inside particle i
                        !-----------------------------------------------
                        !
                        ! P_i = \sum_{c}Fbn_{ij}*a_{ij}/\sum_{c}a_{ij}
                        !
                        !-----------------------------------------------
                        ! total bond area
                        tab(i)  = tab(i) + sb(j, i)                   
                        ! pressure from bonds
                        pb(i)   = pb(i) - p%fbn * sb(j, i)     
                        
                        ! symmetric part
                        tab(j) = tab(j) + sb(j, i)
                        pb(j)  = pb(j) - p%fbn * sb(j, i)
                    end if

                end if

			end if

			! compute sheltering height for particule j on particle i for air and water drag
            ! you have to check both sides of the matrix because it is not symmetric
            if ( shelter .eqv. .true. ) then
                call sheltering(j, i, p)

                ! update local minimum value here because of reduction
                local_hsfa_min_thread(i, thread_id) = min( &
                    local_hsfa_min_thread(i, thread_id), p%hsfa_ji )
                local_hsfa_min_thread(j, thread_id) = min( &
                    local_hsfa_min_thread(j, thread_id), p%hsfa_ij )

                local_hsfw_min_thread(i, thread_id) = min( &
                    local_hsfw_min_thread(i, thread_id), p%hsfw_ji )
                local_hsfw_min_thread(j, thread_id) = min( &
                    local_hsfw_min_thread(j, thread_id), p%hsfw_ij )
            end if

            !-------------------------------------------------------
            !           Computation of diagnotics variables
            !-------------------------------------------------------
            
            if ( flag_diag_stress .eqv. .true. ) then
            ! compute the stress using cauchy stress formula (this needs to be averaged over the size of the particle)
            !-------------------------------------------------------
            !
            !    \sigma_{ij} = 1/A \sum_{c} r_j * Fcn_i
            !
            !-------------------------------------------------------
                block
                    double precision :: force_mag, ri_f, rj_f

                    ! contact forces only count while the pair touches
                    force_mag = sqrt(p%fbn ** 2 + p%fbt ** 2)
                    if ( p%deltan .gt. 0 ) then
                        force_mag = force_mag + &
                                    sqrt(p%fcn ** 2 + p%fct ** 2)
                    end if

                    ri_f = r(i) * force_mag
                    rj_f = r(j) * force_mag

                    sigxx(i) = sigxx(i) - ri_f * p%cosa * p%cosa
                    sigyy(i) = sigyy(i) - ri_f * p%sina * p%sina
                    sigxy(i) = sigxy(i) - ri_f * p%sina * p%cosa
                    sigyx(i) = sigyx(i) - ri_f * p%cosa * p%sina

                    ! Newton's third law equivalent for stress
                    sigxx(j) = sigxx(j) - rj_f * p%cosa * p%cosa
                    sigyy(j) = sigyy(j) - rj_f * p%sina * p%sina
                    sigxy(j) = sigxy(j) - rj_f * p%sina * p%cosa
                    sigyx(j) = sigyx(j) - rj_f * p%cosa * p%sina
                end block
            end if

        end do

        ! verify the bondary conditions for each particle
        if (nx_mask > 0) then
            call verify_bc_mask (i)
        else
            call verify_bc (i)
        end if

    end do
    !$omp end do
    !$omp end parallel

    ! reduce the sheltering coefficient arrays
    if ( shelter .eqv. .true. ) then
        do thread_id = 1, thread_num
            do i = 1, n
                local_hsfa_min(i) = min(local_hsfa_min(i), local_hsfa_min_thread(i, thread_id))
                local_hsfw_min(i) = min(local_hsfw_min(i), local_hsfw_min_thread(i, thread_id))
            end do
        end do

        deallocate(local_hsfa_min_thread)
        deallocate(local_hsfw_min_thread)
    end if

    call timer_off(t_pairs)

    ! reduce the shelter coeff. and the ridged overlap volume
    call timer_on(t_comm)
    call broadcast_shape
    call timer_off(t_comm)

    ! apply the ridging shape changes, identically on every rank
    call timer_on(t_integ)
    call apply_ridging
    call timer_off(t_integ)

    call timer_on(t_forcing)
    !$omp parallel do
    ! compute the total forcing from winds, currents and coriolis
    do i = first_iter, last_iter
        if (.not. active(i)) cycle
        call forcing(i)
        call coriolis(i)
    end do
    !$omp end parallel do
    call timer_off(t_forcing)

    ! reduce all the force variables
    call timer_on(t_comm)
    call force_reduction_fast
    call timer_off(t_comm)

    ! sum all forces together on particule i
    do i = first_iter, last_iter
        tfx_r(i) = fcx_r(i) + fbx_r(i) + fax_r(i) + fwx_r(i) &
                    + fcorx_r(i) + fx_bc_r(i)
        tfy_r(i) = fcy_r(i) + fby_r(i) + fay_r(i) + fwy_r(i) &
                    + fcory_r(i) + fy_bc_r(i)

        ! sum all moments on particule i together
        m_r(i) =  mc_r(i) + mb_r(i) + ma_r(i) + mw_r(i) + m_bc_r(i)

        ! same for stresses
        tsigxx_r(i) = sigxx_r(i) + sigxx_bc_r(i) + sigxx_aw_r(i)
        tsigyy_r(i) = sigyy_r(i) + sigyy_bc_r(i) + sigyy_aw_r(i)
        tsigxy_r(i) = sigxy_r(i) + sigxy_bc_r(i) + sigxy_aw_r(i)
        tsigyx_r(i) = sigyx_r(i) + sigyx_bc_r(i) + sigyx_aw_r(i)

        ! same for pressure
        tp_r(i) =   merge(pc_r(i) / tac_r(i), 0d0, tac_r(i) /= 0d0) &
                  + merge(pb_r(i) / tab_r(i), 0d0, tab_r(i) /= 0d0) &
                  + merge(p_bc_r(i) / ta_bc_r(i), 0d0,              &
                            ta_bc_r(i) /= 0d0)
    end do

    ! broadcast forces to all so that the nodes can each update their x and u
    call timer_on(t_comm)
    call broadcast_total_forces
    call timer_off(t_comm)

    ! deactivate the particles that left the domain, on every rank
    call remove_exited

    ! forces for experiments
!    call normal_forces("ridging", tstep)
!    call gravity

    ! Velocity Verlet: update velocities after force computation
    ! (tstep=1 is Euler initialization; tstep>=2 is Verlet)
    call timer_on(t_integ)
    if ( tstep .ge. 1 ) then
        call velocity
    end if
    call verlet_history
    call timer_off(t_integ)

end subroutine stepper
