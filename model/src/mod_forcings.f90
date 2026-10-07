!=======================================================================
!   Module: wind, current, sheltering, and plate forcings
!
! Allocatable arrays are sized at run time via init_forcings(num_part),
! called from godar.f90 once `n` has been read from the input file.
!=======================================================================

module forcings

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! ambient (uniform) wind and currents
    real(dp) ::                 &
        ua       ,              & ! base wind in x                [m/s]
        va       ,              & ! base wind in y                [m/s]
        uw       ,              & ! base currents in x            [m/s]
        vw                        ! base currents in y            [m/s]

    real(dp), allocatable ::    &
        ua_i (:),               & ! per-particle wind in x        [m/s]
        va_i (:),               & ! per-particle wind in y        [m/s]
        uw_i (:),               & ! per-particle currents in x    [m/s]
        vw_i (:)                  ! per-particle currents in y    [m/s]

    ! forcing forces and torques
    real(dp), allocatable ::    &
        fax (:),                & ! wind force in x                 [N]
        fay (:),                & ! wind force in y                 [N]
        fwx (:),                & ! water force in x                [N]
        fwy (:),                & ! water force in y                [N]
        mw  (:),                & ! water drag moment             [N*m]
        ma  (:)                   ! air drag moment               [N*m]

    ! plate boundary forces (legacy rectangular path)
    real(dp) ::                 &
        pfn,                    & ! normal force on the plates
        pfs                       ! shear force on the plates

contains

    subroutine alloc_forcings(num_part)

        integer, intent(in) :: num_part

        ! zero-initialised, as the old common blocks were

        ! per-particle ambient wind and currents
        allocate(ua_i(num_part), source=0.0_dp)
        allocate(va_i(num_part), source=0.0_dp)
        allocate(uw_i(num_part), source=0.0_dp)
        allocate(vw_i(num_part), source=0.0_dp)

        ! forcing forces and torques
        allocate(fax(num_part), source=0.0_dp)
        allocate(fay(num_part), source=0.0_dp)
        allocate(fwx(num_part), source=0.0_dp)
        allocate(fwy(num_part), source=0.0_dp)
        allocate(mw (num_part), source=0.0_dp)
        allocate(ma (num_part), source=0.0_dp)

    end subroutine alloc_forcings


    subroutine dealloc_forcings

        ! per-particle ambient wind and currents
        deallocate(ua_i)
        deallocate(va_i)
        deallocate(uw_i)
        deallocate(vw_i)

        ! forcing forces and torques
        deallocate(fax)
        deallocate(fay)
        deallocate(fwx)
        deallocate(fwy)
        deallocate(mw )
        deallocate(ma )

    end subroutine dealloc_forcings

end module forcings
