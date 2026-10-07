!=======================================================================
!   Module: dynamic particle state, forces, and pair decomposition
!
! Allocatable arrays are sized at run time via init_variables(num_part),
! called from godar.f90 once `n` has been read from the input file.
!=======================================================================

module variables

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! particle state
    real(dp), allocatable ::    &
        x      (:),             & ! x position                      [m]
        y      (:),             & ! y position                      [m]
        r      (:),             & ! radius                          [m]
        h      (:),             & ! thickness                       [m]
        hfa    (:),             & ! freeboard                       [m]
        hfw    (:),             & ! draft                           [m]
        mass   (:),             & ! mass                           [kg]
        inertia(:),             & ! moment of inertia          [kg*m^2]
        u      (:),             & ! velocity in x                 [m/s]
        v      (:),             & ! velocity in y                 [m/s]
        theta  (:),             & ! angular position              [rad]
        omega  (:)                ! angular velocity            [rad/s]

    ! force fields
    real(dp), allocatable ::    &
        fcx  (:),               & ! summed contact force in x       [N]
        fcy  (:),               & ! summed contact force in y       [N]
        mc   (:),               & ! moment due to contact         [N*m]
        tfx  (:),               & ! total force in x                [N]
        tfy  (:),               & ! total force in y                [N]
        m    (:),               & ! total moment                  [N*m]
        fcorx(:),               & ! Coriolis force in x             [N]
        fcory(:)                  ! Coriolis force in y             [N]

    ! contact history per pair (j,i), kept between steps; the values
    ! computed for a pair during a step are in pair_t (mod_pair.f90)
    real(dp), allocatable ::    &
        deltat   (:,:),         & ! accumulated tangent overlap     [m]
        thetarelc(:,:)            ! relative angular position     [rad]

    ! boundary-condition forces and per-particle BC history
    real(dp), allocatable ::    &
        fn_bc       (:),        & ! normal force                  [N]
        ft_bc       (:),        & ! tangential force              [N]
        fr_bc       (:),        & ! rolling friction force        [N]
        fx_bc       (:),        & ! BC force in x                 [N]
        fy_bc       (:),        & ! BC force in y                 [N]
        mc_bc       (:),        & ! BC rolling moment           [N*m]
        m_bc        (:),        & ! total BC moment             [N*m]
        deltat_bc1  (:),        & ! tangential compression, bd 1  [m]
        deltat_bc2  (:),        & ! tangential compression, bd 2  [m]
        theta_bc1   (:),        & ! rolling angle, bd 1         [rad]
        theta_bc2   (:),        & ! rolling angle, bd 2         [rad]
        delt_ridge_bc(:)          ! tangent overlap for BC ridging[m]

    ! Velocity-Verlet stored acceleration history
    real(dp), allocatable ::    &
        ax_nm1      (:),        & ! x acceleration at n-1     [m/s^2]
        ay_nm1      (:),        & ! y acceleration at n-1     [m/s^2]
        atheta_nm1  (:)           ! angular accel at n-1    [rad/s^2]

    ! activity flag: .true. while particle is inside the (mask) grid
    logical, allocatable :: active(:)

    ! exit flag: set when a particle leaves the domain during a step,
    ! applied to active after the force loop (remove_exited)
    logical, allocatable :: exited(:)

    ! overlap volume converted into thickness after the force loop
    ! (apply_ridging), mass conserving                          [m^3]
    real(dp), allocatable :: dvol(:)

contains

    subroutine alloc_variables(num_part)

        integer, intent(in) :: num_part

        ! zero-initialised, as the old common blocks were

        ! particle state
        allocate(x      (num_part), source=0.0_dp)
        allocate(y      (num_part), source=0.0_dp)
        allocate(r      (num_part), source=0.0_dp)
        allocate(h      (num_part), source=0.0_dp)
        allocate(hfa    (num_part), source=0.0_dp)
        allocate(hfw    (num_part), source=0.0_dp)
        allocate(mass   (num_part), source=0.0_dp)
        allocate(inertia(num_part), source=0.0_dp)
        allocate(u      (num_part), source=0.0_dp)
        allocate(v      (num_part), source=0.0_dp)
        allocate(theta  (num_part), source=0.0_dp)
        allocate(omega  (num_part), source=0.0_dp)

        ! force fields
        allocate(fcx  (num_part), source=0.0_dp)
        allocate(fcy  (num_part), source=0.0_dp)
        allocate(mc   (num_part), source=0.0_dp)
        allocate(tfx  (num_part), source=0.0_dp)
        allocate(tfy  (num_part), source=0.0_dp)
        allocate(m    (num_part), source=0.0_dp)
        allocate(fcorx(num_part), source=0.0_dp)
        allocate(fcory(num_part), source=0.0_dp)

        ! contact history per pair
        allocate(deltat   (num_part, num_part), source=0.0_dp)
        allocate(thetarelc(num_part, num_part), source=0.0_dp)

        ! boundary-condition forces and BC history
        allocate(fn_bc        (num_part), source=0.0_dp)
        allocate(ft_bc        (num_part), source=0.0_dp)
        allocate(fr_bc        (num_part), source=0.0_dp)
        allocate(fx_bc        (num_part), source=0.0_dp)
        allocate(fy_bc        (num_part), source=0.0_dp)
        allocate(mc_bc        (num_part), source=0.0_dp)
        allocate(m_bc         (num_part), source=0.0_dp)
        allocate(deltat_bc1   (num_part), source=0.0_dp)
        allocate(deltat_bc2   (num_part), source=0.0_dp)
        allocate(theta_bc1    (num_part), source=0.0_dp)
        allocate(theta_bc2    (num_part), source=0.0_dp)
        allocate(delt_ridge_bc(num_part), source=0.0_dp)

        ! Velocity-Verlet stored acceleration history
        allocate(ax_nm1    (num_part), source=0.0_dp)
        allocate(ay_nm1    (num_part), source=0.0_dp)
        allocate(atheta_nm1(num_part), source=0.0_dp)

        ! activity flag
        allocate(active(num_part), source=.false.)

        ! exit flag
        allocate(exited(num_part), source=.false.)

        ! ridged overlap volume
        allocate(dvol(num_part), source=0.0_dp)

    end subroutine alloc_variables


    subroutine dealloc_variables

        ! particle state
        deallocate(x      )
        deallocate(y      )
        deallocate(r      )
        deallocate(h      )
        deallocate(hfa    )
        deallocate(hfw    )
        deallocate(mass   )
        deallocate(inertia)
        deallocate(u      )
        deallocate(v      )
        deallocate(theta  )
        deallocate(omega  )

        ! force fields
        deallocate(fcx  )
        deallocate(fcy  )
        deallocate(mc   )
        deallocate(tfx  )
        deallocate(tfy  )
        deallocate(m    )
        deallocate(fcorx)
        deallocate(fcory)

        ! contact history per pair
        deallocate(deltat   )
        deallocate(thetarelc)

        ! boundary-condition forces and BC history
        deallocate(fn_bc        )
        deallocate(ft_bc        )
        deallocate(fr_bc        )
        deallocate(fx_bc        )
        deallocate(fy_bc        )
        deallocate(mc_bc        )
        deallocate(m_bc         )
        deallocate(deltat_bc1   )
        deallocate(deltat_bc2   )
        deallocate(theta_bc1    )
        deallocate(theta_bc2    )
        deallocate(delt_ridge_bc)

        ! Velocity-Verlet stored acceleration history
        deallocate(ax_nm1    )
        deallocate(ay_nm1    )
        deallocate(atheta_nm1)

        ! activity flag
        deallocate(active)

        ! exit flag
        deallocate(exited)

        ! ridged overlap volume
        deallocate(dvol)

    end subroutine dealloc_variables

end module variables
