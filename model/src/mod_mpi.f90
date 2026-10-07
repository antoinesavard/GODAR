!=======================================================================
!   Module: MPI rank state and per-rank reduction buffers
!
! Allocatable buffers are sized at run time via init_mpi_buffers(num_part),
! called from godar.f90 once `n` has been read from the input file.
! Module storage (rather than the old CB_mpi.h common block) also avoids
! the GOTPCREL relocation issue seen with large `n` on some toolchains.
!=======================================================================

module mpi_var

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! rank-level state
    integer ::      &
        rank      , & ! this process' rank
        n_ranks   , & ! total number of processes
        ierr          ! generic MPI error code

    ! local-slice info on this rank
    integer ::      &
        local_n   , & ! number of particles owned by this rank
        local_disp, & ! displacement (zero-based) in the global vector
        first_iter, & ! first particle index this rank iterates
        last_iter     ! last  particle index this rank iterates

    ! reduction receive buffers (per-particle aggregates after reduce)
    real(dp), allocatable ::    &
        ! contact
        fcx_r       (:),        &
        fcy_r       (:),        &
        mc_r        (:),        &
        ! bond
        fbx_r       (:),        &
        fby_r       (:),        &
        mb_r        (:),        &
        ! forcing
        fwx_r       (:),        &
        fwy_r       (:),        &
        mw_r        (:),        &
        fax_r       (:),        &
        fay_r       (:),        &
        ma_r        (:),        &
        fcorx_r     (:),        &
        fcory_r     (:),        &
        ! boundary
        fx_bc_r     (:),        &
        fy_bc_r     (:),        &
        m_bc_r      (:),        &
        ! stress (pair)
        sigxx_r     (:),        &
        sigyy_r     (:),        &
        sigxy_r     (:),        &
        sigyx_r     (:),        &
        ! stress (boundary)
        sigxx_bc_r  (:),        &
        sigyy_bc_r  (:),        &
        sigxy_bc_r  (:),        &
        sigyx_bc_r  (:),        &
        ! stress (atm + water)
        sigxx_aw_r  (:),        &
        sigyy_aw_r  (:),        &
        sigxy_aw_r  (:),        &
        sigyx_aw_r  (:),        &
        ! pressure
        tac_r       (:),        &
        tab_r       (:),        &
        pc_r        (:),        &
        pb_r        (:),        &
        ! pressure (boundary)
        ta_bc_r     (:),        &
        p_bc_r      (:)

    ! global totals (filled after broadcast_total_forces)
    real(dp), allocatable ::    &
        m_r         (:),        &
        tfx_r       (:),        &
        tfy_r       (:),        &
        tsigxx_r    (:),        &
        tsigyy_r    (:),        &
        tsigxy_r    (:),        &
        tsigyx_r    (:),        &
        tp_r        (:)

    ! sheltering reduction buffers
    real(dp), allocatable ::    &
        hsfa_min_r    (:),      & ! global min of hsfa per particle
        hsfw_min_r    (:),      & ! global min of hsfw per particle
        local_hsfa_min(:),      & ! local min on this rank
        local_hsfw_min(:)         ! local min on this rank

contains

    subroutine alloc_mpi(num_part)

        integer, intent(in) :: num_part

        ! zero-initialised, as the old common blocks were

        ! reduction receive buffers (contact)
        allocate(fcx_r(num_part), source=0.0_dp)
        allocate(fcy_r(num_part), source=0.0_dp)
        allocate(mc_r (num_part), source=0.0_dp)

        ! reduction receive buffers (bond)
        allocate(fbx_r(num_part), source=0.0_dp)
        allocate(fby_r(num_part), source=0.0_dp)
        allocate(mb_r (num_part), source=0.0_dp)

        ! reduction receive buffers (forcing)
        allocate(fwx_r  (num_part), source=0.0_dp)
        allocate(fwy_r  (num_part), source=0.0_dp)
        allocate(mw_r   (num_part), source=0.0_dp)
        allocate(fax_r  (num_part), source=0.0_dp)
        allocate(fay_r  (num_part), source=0.0_dp)
        allocate(ma_r   (num_part), source=0.0_dp)
        allocate(fcorx_r(num_part), source=0.0_dp)
        allocate(fcory_r(num_part), source=0.0_dp)

        ! reduction receive buffers (boundary)
        allocate(fx_bc_r(num_part), source=0.0_dp)
        allocate(fy_bc_r(num_part), source=0.0_dp)
        allocate(m_bc_r (num_part), source=0.0_dp)

        ! reduction receive buffers (stress, pair)
        allocate(sigxx_r(num_part), source=0.0_dp)
        allocate(sigyy_r(num_part), source=0.0_dp)
        allocate(sigxy_r(num_part), source=0.0_dp)
        allocate(sigyx_r(num_part), source=0.0_dp)

        ! reduction receive buffers (stress, boundary)
        allocate(sigxx_bc_r(num_part), source=0.0_dp)
        allocate(sigyy_bc_r(num_part), source=0.0_dp)
        allocate(sigxy_bc_r(num_part), source=0.0_dp)
        allocate(sigyx_bc_r(num_part), source=0.0_dp)

        ! reduction receive buffers (stress, atm + water)
        allocate(sigxx_aw_r(num_part), source=0.0_dp)
        allocate(sigyy_aw_r(num_part), source=0.0_dp)
        allocate(sigxy_aw_r(num_part), source=0.0_dp)
        allocate(sigyx_aw_r(num_part), source=0.0_dp)

        ! reduction receive buffers (pressure)
        allocate(tac_r  (num_part), source=0.0_dp)
        allocate(tab_r  (num_part), source=0.0_dp)
        allocate(pc_r   (num_part), source=0.0_dp)
        allocate(pb_r   (num_part), source=0.0_dp)
        allocate(ta_bc_r(num_part), source=0.0_dp)
        allocate(p_bc_r (num_part), source=0.0_dp)

        ! global totals
        allocate(m_r     (num_part), source=0.0_dp)
        allocate(tfx_r   (num_part), source=0.0_dp)
        allocate(tfy_r   (num_part), source=0.0_dp)
        allocate(tsigxx_r(num_part), source=0.0_dp)
        allocate(tsigyy_r(num_part), source=0.0_dp)
        allocate(tsigxy_r(num_part), source=0.0_dp)
        allocate(tsigyx_r(num_part), source=0.0_dp)
        allocate(tp_r    (num_part), source=0.0_dp)

        ! sheltering reduction buffers
        allocate(hsfa_min_r    (num_part), source=0.0_dp)
        allocate(hsfw_min_r    (num_part), source=0.0_dp)
        allocate(local_hsfa_min(num_part), source=0.0_dp)
        allocate(local_hsfw_min(num_part), source=0.0_dp)

    end subroutine alloc_mpi


    subroutine dealloc_mpi

        ! reduction receive buffers (contact)
        deallocate(fcx_r)
        deallocate(fcy_r)
        deallocate(mc_r )

        ! reduction receive buffers (bond)
        deallocate(fbx_r)
        deallocate(fby_r)
        deallocate(mb_r )

        ! reduction receive buffers (forcing)
        deallocate(fwx_r  )
        deallocate(fwy_r  )
        deallocate(mw_r   )
        deallocate(fax_r  )
        deallocate(fay_r  )
        deallocate(ma_r   )
        deallocate(fcorx_r)
        deallocate(fcory_r)

        ! reduction receive buffers (boundary)
        deallocate(fx_bc_r)
        deallocate(fy_bc_r)
        deallocate(m_bc_r )

        ! reduction receive buffers (stress, pair)
        deallocate(sigxx_r)
        deallocate(sigyy_r)
        deallocate(sigxy_r)
        deallocate(sigyx_r)

        ! reduction receive buffers (stress, boundary)
        deallocate(sigxx_bc_r)
        deallocate(sigyy_bc_r)
        deallocate(sigxy_bc_r)
        deallocate(sigyx_bc_r)

        ! reduction receive buffers (stress, atm + water)
        deallocate(sigxx_aw_r)
        deallocate(sigyy_aw_r)
        deallocate(sigxy_aw_r)
        deallocate(sigyx_aw_r)

        ! reduction receive buffers (pressure)
        deallocate(tac_r  )
        deallocate(tab_r  )
        deallocate(pc_r   )
        deallocate(pb_r   )
        deallocate(ta_bc_r)
        deallocate(p_bc_r )

        ! global totals
        deallocate(m_r     )
        deallocate(tfx_r   )
        deallocate(tfy_r   )
        deallocate(tsigxx_r)
        deallocate(tsigyy_r)
        deallocate(tsigxy_r)
        deallocate(tsigyx_r)
        deallocate(tp_r    )

        ! sheltering reduction buffers
        deallocate(hsfa_min_r    )
        deallocate(hsfw_min_r    )
        deallocate(local_hsfa_min)
        deallocate(local_hsfw_min)

    end subroutine dealloc_mpi

end module mpi_var
