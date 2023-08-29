subroutine stepper (tstep)

    use omp_lib
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_KdTree, only: KdTree, KdTreeSearch
    use dArgDynamicArray_Class, only: dArgDynamicArray
    use m_strings, only: str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"
    include "CB_forcings.h"

    integer :: i, j, k
    integer, intent(in) :: tstep
    type(KdTree) :: tree
    type(KdTreeSearch) :: search
    type(dArgDynamicArray) :: da

    ! Build the tree
    tree = KdTree(x, y)
        
    ! reinitialize force arrays for contact and bonds
    do i = 1, n
        mc(i)    = 0d0
        mb(i)    = 0d0
        fcx(i)   = 0d0
        fcy(i)   = 0d0
        fbx(i)   = 0d0
        fby(i)   = 0d0
        ! and total force arrays
        m(i)   = 0d0
		tfx(i) = 0d0
        tfy(i) = 0d0
    end do
    
    ! put yourself in the referential of the ith particle
	! loop through all j particles and compute interactions

    !$omp parallel do schedule(guided) &
    !$omp private(i,j,da) &
    !$omp reduction(+:fcx,fcy,fbx,fby,mc,mb)
    do i = n, 1, -1
        ! Find all the particles j near i
        da = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), &
                            radius = r(i) + rtree)
        ! loop over the nearest neighbors except the first because this is the particle i
        do k = 1, size(da%i%values) - 1
            j = da%i%values(k + 1)

            ! only compute lower triangular matrices
            if (i .ge. j) then
                cycle
            end if

			! compute relative position and velocity
            call rel_pos_vel (j, i)

			! bond initialization
			if ( tstep .eq. 1 ) then
				if ( -deltan(j, i) .le. 5d-1 * r(i)) then ! can be fancier
					bond (j, i) = 1
                end if
                call bond_properties (j, i)
			end if

            ! verify if two particles are colliding
            if ( deltan(j,i) .gt. 0 ) then

                call contact_forces (j, i)
				call bond_creation (j, i) ! to implement

				! update contact force on particle i by particle j
                fcx(i) = fcx(i) - fcn(j,i) * cosa(j,i)
                fcy(i) = fcy(i) - fcn(j,i) * sina(j,i)

                ! update moment on particule i by particule j due to tangent contact 
                mc(i) = mc(i) - r(i) * fct(j,i) - mcc(j,i)

                ! Newton's third law
                ! update contact force on particle j by particle i
                fcx(j) = fcx(j) + fcn(j,i) * cosa(j,i)
                fcy(j) = fcy(j) + fcn(j,i) * sina(j,i)

                ! update moment on particule j by particule i due to tangent contact 
                mc(j) = mc(j) - r(j) * fct(j,i) + mcc(j,i)

            end if

			! compute forces from bonds between particle i and j
			if ( bond (j, i) .eq. 1 ) then

				call bond_forces (j, i)
				call bond_breaking (j, i)

                if ( bond (j, i) .eq. 1 ) then
                    ! update force on particle i by j due to bond
                    fbx(i) = fbx(i) - fbn(j,i) * cosa(j,i) +    &
                                        fbt(j,i) * sina(j,i)
                    fby(i) = fby(i) - fbn(j,i) * sina(j,i) -    &
                                        fbt(j,i) * cosa(j,i)

                    ! update moment on particule i by j to to bond
                    mb(i) = mb(i) - r(i) * fbt(j,i) - mbb(j, i)

                    ! Newton's third law
                    ! update force on particle j by i due to bond
                    fbx(j) = fbx(j) + fbn(j,i) * cosa(j,i) -    &
                                        fbt(j,i) * sina(j,i)
                    fby(j) = fby(j) + fbn(j,i) * sina(j,i) +    &
                                        fbt(j,i) * cosa(j,i)
    
                    ! update moment on particule j by i due to bond
                    mb(j) = mb(j) - r(j) * fbt(j,i) + mbb(j, i)
                end if

			end if

			! compute sheltering height for particule j on particle i for air and water drag
            call sheltering(j, i)

        end do

        ! compute the total forcing from winds, currents and coriolis on particule i
        call forcing (i)
        call coriolis(i)

    end do
    !$omp end parallel do

    ! deallocate tree memory
    call tree%deallocate()

    ! sum all forces together on particule i
    do i = 1, n
        tfx(i) = fcx(i) + fbx(i) + fax(i) + fwx(i) + fcorx(i)
        tfy(i) = fcy(i) + fby(i) + fay(i) + fwy(i) + fcory(i)

        ! sum all moments on particule i together
        m(i) =  mc(i) + mb(i) + ma(i) + mw(i)
    end do

    ! forces on side particles for experiments
    !call experiment_forces

	! integration in time
    call velocity
    call euler

end subroutine stepper


subroutine experiment_forces

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"
    include "CB_forcings.h"


    tfx(44) = 0d0
    tfy(44) = 0d0
    tfx(45) = 0d0
    tfy(45) = 0d0

    ! bottom plate on the left
    !tfy(1) = tfy(1) + 1d8
    !tfy(2) = tfy(2) + 1d8
    !tfy(3) = tfy(3) + 1d8

    ! top plate on the right
    tfy(6) = tfy(6) - 5d7
    tfy(4) = tfy(4) - 5d7
    tfy(7) = tfy(7) - 5d7

    ! left plate
    !tfx(2) = tfx(2) + 1d8
    !tfx(8) = tfx(8) + 1d8
    !tfx(14) = tfx(14) + 1d8
    !tfx(20) = tfx(20) + 1d8
    !tfx(26) = tfx(26) + 1d8

    ! right plate
    tfx(7) = tfx(7) - 5d7
    tfx(13) = tfx(13) - 5d7
    tfx(19) = tfx(19) - 5d7
    tfx(25) = tfx(25) - 5d7
    tfx(31) = tfx(31) - 5d7

end subroutine experiment_forces
