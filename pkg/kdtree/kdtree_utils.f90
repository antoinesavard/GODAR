module kdtree_utils

    use global_KdTree
    use m_KdTree, only: KdTree, KdTreeSearch
    use dArgDynamicArray_Class, only: dArgDynamicArray

    implicit none

    ! neighbours j > i of each particle i: the particles within its
    ! search radius when the lists were built, by increasing distance
    type :: neighbour_list
        integer, allocatable :: j(:)
    end type neighbour_list

    type(neighbour_list), allocatable :: neighbours(:)

    ! positions of the particles when the lists were built
    double precision, allocatable :: x_built(:), y_built(:)

contains

    subroutine tree_building(tstep, ntree, xtree, ytree)

        integer, intent(in) :: tstep, ntree
        double precision, intent(in) :: xtree(:), ytree(:)

        ! Build the tree
        ! you can skip any timesteps using a copied array of the positions as the tree building process does not make a copy of the data and changing the data on which the tree is based will throw a seg fault.
        if ( mod(tstep, ntree) == 0 ) then
            tree = KdTree(xtree, ytree)
            is_built = .true.
        end if

    end subroutine tree_building


    subroutine tree_cleanup(tstep, ntree)

        integer, intent(in) :: tstep, ntree

        if ( mod(tstep, ntree) == ntree - 1 ) then
            call tree%deallocate()
            is_built = .false.
        end if

    end subroutine tree_cleanup


    subroutine build_neighbours(x, y, radius, first, last)

        ! builds the neighbour lists of the particles first to last with
        ! a kd-tree of the current positions, then frees the tree

        double precision, intent(in) :: x(:), y(:), radius(:)
        integer, intent(in) :: first, last

        type(KdTree) :: tree_nb
        type(KdTreeSearch) :: search
        type(dArgDynamicArray) :: da
        integer :: i

        if (.not. allocated(neighbours)) allocate(neighbours(size(x)))

        tree_nb = KdTree(x, y)

        !$omp parallel do private(da) schedule(dynamic, 16)
        do i = first, last
            da = search%kNearest(tree_nb, x, y, xQuery = x(i), &
                                 yQuery = y(i), radius = radius(i))
            neighbours(i)%j = pack(da%i%values, da%i%values > i)
        end do
        !$omp end parallel do

        call tree_nb%deallocate()

        x_built = x
        y_built = y

    end subroutine build_neighbours


    logical function neighbours_outdated(x, y, margin)

        ! true when a particle moved more than margin / 2 since the lists
        ! were built: two particles may then have come closer by more than
        ! margin, so a pair that is not listed could be in contact

        double precision, intent(in) :: x(:), y(:), margin

        if ( .not. allocated(x_built) .or. margin .le. 0d0 ) then
            neighbours_outdated = .true.
        else
            neighbours_outdated = &
                4 * maxval( (x - x_built)**2 + (y - y_built)**2 ) > margin**2
        end if

    end function neighbours_outdated

end module kdtree_utils


! subroutine kdtree_update()

!     implicit none

!     include "parameter.h"
!     include "CB_variables.h"
!     include "CB_const.h"

!     double precision :: max_veln, max_velt, max_veln_bc, max_velt_bc
!     double precision :: min_r1, min_r2

!     max_veln = maxval( veln )
!     max_velt = maxval( velt )
!     max_veln_bc = maxval( veln_bc )
!     max_velt_bc = maxval( velt_bc )

!     max_vel = max( max_veln, max_velt, max_veln_bc, max_velt_bc)

!     min_r1 = minval( r )
!     min_r2 = minval( r, mask = r > min_r1)

! end subroutine kdtree_update
