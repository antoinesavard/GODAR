!=======================================================================
!   Module: orchestrate allocation / deallocation of all per-particle
!           module arrays
!
! Each data module exposes its own alloc_<module>(num_part) and
! dealloc_<module>() routines. This module just chains them together so
! godar.f90 can do `call allocate_all(n)` after reading n from input.
!=======================================================================

module alloc

    use parameters,  only: n
    use variables,   only: alloc_variables,   dealloc_variables
    use bonds,       only: alloc_bonds,       dealloc_bonds
    use diagnostics, only: alloc_diagnostics, dealloc_diagnostics
    use forcings,    only: alloc_forcings,    dealloc_forcings
    use mpi_var,     only: alloc_mpi,         dealloc_mpi

    implicit none

contains

    subroutine allocate_all(num_part)

        integer, intent(in) :: num_part

        call alloc_variables  (num_part)
        call alloc_bonds      (num_part)
        call alloc_diagnostics(num_part)
        call alloc_forcings   (num_part)
        call alloc_mpi        (num_part)

    end subroutine allocate_all


    subroutine deallocate_all

        call dealloc_variables
        call dealloc_bonds
        call dealloc_diagnostics
        call dealloc_forcings
        call dealloc_mpi

    end subroutine deallocate_all

end module alloc
