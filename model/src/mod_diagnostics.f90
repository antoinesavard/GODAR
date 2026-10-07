!=======================================================================
!   Module: diagnostic variables (stress, pressure, areas)
!
! Allocatable arrays are sized at run time via init_diagnostics(num_part),
! called from godar.f90 once `n` has been read from the input file.
!=======================================================================

module diagnostics

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! stress accumulators
    real(dp), allocatable ::    &
        sigxx    (:),           & ! pair-contact + bond contributions
        sigyy    (:),           &
        sigxy    (:),           &
        sigyx    (:),           &
        ! totals
        tsigxx   (:),           &
        tsigyy   (:),           &
        tsigxy   (:),           &
        tsigyx   (:),           &
        ! boundary contributions
        sigxx_bc (:),           &
        sigyy_bc (:),           &
        sigxy_bc (:),           &
        sigyx_bc (:),           &
        ! atmosphere + water contributions
        sigxx_aw (:),           &
        sigyy_aw (:),           &
        sigxy_aw (:),           &
        sigyx_aw (:)

    ! pressure/contact-area accumulators
    real(dp), allocatable ::    &
        tac   (:),              & ! total area of contacts
        tab   (:),              & ! total area of bonds
        pc    (:),              & ! pressure from contacts
        pb    (:),              & ! pressure from bonds
        tp    (:),              & ! total pressure
        ! boundary contributions
        ta_bc (:),              & ! total area of BC contacts
        p_bc  (:)                 ! total pressure from BC contacts

contains

    subroutine alloc_diagnostics(num_part)

        integer, intent(in) :: num_part

        ! zero-initialised, as the old common blocks were

        ! stress accumulators
        allocate(sigxx   (num_part), source=0.0_dp)
        allocate(sigyy   (num_part), source=0.0_dp)
        allocate(sigxy   (num_part), source=0.0_dp)
        allocate(sigyx   (num_part), source=0.0_dp)
        allocate(tsigxx  (num_part), source=0.0_dp)
        allocate(tsigyy  (num_part), source=0.0_dp)
        allocate(tsigxy  (num_part), source=0.0_dp)
        allocate(tsigyx  (num_part), source=0.0_dp)
        allocate(sigxx_bc(num_part), source=0.0_dp)
        allocate(sigyy_bc(num_part), source=0.0_dp)
        allocate(sigxy_bc(num_part), source=0.0_dp)
        allocate(sigyx_bc(num_part), source=0.0_dp)
        allocate(sigxx_aw(num_part), source=0.0_dp)
        allocate(sigyy_aw(num_part), source=0.0_dp)
        allocate(sigxy_aw(num_part), source=0.0_dp)
        allocate(sigyx_aw(num_part), source=0.0_dp)

        ! pressure/contact-area accumulators
        allocate(tac  (num_part), source=0.0_dp)
        allocate(tab  (num_part), source=0.0_dp)
        allocate(pc   (num_part), source=0.0_dp)
        allocate(pb   (num_part), source=0.0_dp)
        allocate(tp   (num_part), source=0.0_dp)
        allocate(ta_bc(num_part), source=0.0_dp)
        allocate(p_bc (num_part), source=0.0_dp)

    end subroutine alloc_diagnostics


    subroutine dealloc_diagnostics

        ! stress accumulators
        deallocate(sigxx   )
        deallocate(sigyy   )
        deallocate(sigxy   )
        deallocate(sigyx   )
        deallocate(tsigxx  )
        deallocate(tsigyy  )
        deallocate(tsigxy  )
        deallocate(tsigyx  )
        deallocate(sigxx_bc)
        deallocate(sigyy_bc)
        deallocate(sigxy_bc)
        deallocate(sigyx_bc)
        deallocate(sigxx_aw)
        deallocate(sigyy_aw)
        deallocate(sigxy_aw)
        deallocate(sigyx_aw)

        ! pressure/contact-area accumulators
        deallocate(tac  )
        deallocate(tab  )
        deallocate(pc   )
        deallocate(pb   )
        deallocate(tp   )
        deallocate(ta_bc)
        deallocate(p_bc )

    end subroutine dealloc_diagnostics

end module diagnostics
