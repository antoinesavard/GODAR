subroutine tree_building(tstep, xtree, ytree)

    use global_KdTree
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_KdTree, only: KdTree
    use dArgDynamicArray_Class, only: dArgDynamicArray
    use m_strings, only: str
    
    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: tstep
    double precision, dimension(n), intent(in) :: xtree, ytree

    ! Build the tree
    ! you can skip any timesteps using a copied array of the positions as the tree building process does not make a copy of the data and changing the data on which the tree is based will throw a seg fault.
    if ( tstep == 1 .or. mod(tstep, int(ntree)) == 0 ) then
        tree = KdTree(xtree, ytree)
        is_built = .true.
    end if


end subroutine tree_building


subroutine tree_cleanup(tstep)

    use global_KdTree
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_KdTree, only: KdTree
    use dArgDynamicArray_Class, only: dArgDynamicArray
    use m_strings, only: str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: tstep

    if ( mod(tstep, int(ntree)) == ntree - 1 ) then
        call tree%deallocate()
        is_built = .false.
    end if


end subroutine tree_cleanup


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