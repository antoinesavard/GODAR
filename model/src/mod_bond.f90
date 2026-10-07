!=======================================================================
!   Module: bond data (state, forces, stresses, damage)
!
! Allocatable arrays are sized at run time via init_bonds(num_part),
! called from godar.f90 once `n` has been read from the input file.
!=======================================================================

module bonds

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! bond diagnostics
    integer ::                  &
        num_bonds        ,      &
        num_bonds_nm1    ,      &
        num_bonds_broken

    ! global physical parameters
    real(dp) ::         &
        eb           ,  & ! elastic stiffness
        lambda_rb    ,  & ! radius parameter
        lambda_lb    ,  & ! length parameter
        sigmacb_crit ,  & ! critical compressive strength
        tau_crit     ,  & ! critical shear strength
        bond_lim     ,  & ! distance limit between particles for bond
        dmax         ,  & ! max damage for bond breaking
        dtd          ,  & ! damage evolution timescale          [s]
        dth          ,  & ! healing timescale                   [s]
        phi_int           ! internal friction angle           [deg]

    ! individual physical properties, kept between steps; the forces
    ! and deformation of a bond during a step are in pair_t
    real(dp), allocatable ::    &
        rb      (:,:),          & ! bond half width
        hb      (:,:),          & ! bond thickness
        lb      (:,:),          & ! bond length
        sb      (:,:),          & ! bond cross-sectional area
        ib      (:,:),          & ! bond moment of inertia
        damageb (:,:)             ! damage variable for bond i-j

    ! angle variables for bending
    real(dp), allocatable ::    &
        theta_offset(:,:),      & ! angle offset from lab frame
        cosa_old    (:,:),      & ! old value of cosa
        sina_old    (:,:),      & ! old value of sina
        alpha_total (:,:)         ! total beam angle

    ! summed forces in bonds
    real(dp), allocatable ::    &
        mb  (:),                & ! total moment due to bonds on i
        fbx (:),                & ! total bond force in x on i
        fby (:)                   ! total bond force in y on i

    ! bond presence (0 or 1)
    integer, allocatable :: bond(:,:)

contains

    subroutine alloc_bonds(num_part)

        integer, intent(in) :: num_part

        ! zero-initialised, as the old common blocks were

        ! individual physical properties
        allocate(rb       (num_part, num_part), source=0.0_dp)
        allocate(hb       (num_part, num_part), source=0.0_dp)
        allocate(lb       (num_part, num_part), source=0.0_dp)
        allocate(sb       (num_part, num_part), source=0.0_dp)
        allocate(ib       (num_part, num_part), source=0.0_dp)
        allocate(damageb  (num_part, num_part), source=0.0_dp)

        ! angle variables for bending
        allocate(theta_offset(num_part, num_part), source=0.0_dp)
        allocate(cosa_old   (num_part, num_part), source=0.0_dp)
        allocate(sina_old   (num_part, num_part), source=0.0_dp)
        allocate(alpha_total(num_part, num_part), source=0.0_dp)

        ! summed forces in bonds
        allocate(mb (num_part), source=0.0_dp)
        allocate(fbx(num_part), source=0.0_dp)
        allocate(fby(num_part), source=0.0_dp)

        ! bond presence
        allocate(bond(num_part, num_part), source=0)

    end subroutine alloc_bonds


    subroutine dealloc_bonds

        ! individual physical properties
        deallocate(rb       )
        deallocate(hb       )
        deallocate(lb       )
        deallocate(sb       )
        deallocate(ib       )
        deallocate(damageb  )

        ! angle variables for bending
        deallocate(theta_offset)
        deallocate(cosa_old   )
        deallocate(sina_old   )
        deallocate(alpha_total)

        ! summed forces in bonds
        deallocate(mb )
        deallocate(fbx)
        deallocate(fby)

        ! bond presence
        deallocate(bond)

    end subroutine dealloc_bonds

end module bonds
