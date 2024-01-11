subroutine stepper (tstep)

    use mpi
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
    include "CB_options.h"
    include "CB_mpi.h"

    integer :: i, j, k
    integer, intent(in) :: tstep
    type(KdTree) :: tree
    type(KdTreeSearch) :: search
    type(dArgDynamicArray) :: da

    ! Build the tree
    tree = KdTree(x, y)

    call reset_forces
    
    ! put yourself in the referential of the ith particle
	! loop through all j particles and compute interactions

    !$omp parallel do schedule(dynamic, 1) &
    !$omp private(i,j,da) &
    !$omp reduction(+:fcx,fcy,fbx,fby,mc,mb)
    do i = last_iter, first_iter, -1
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
			if ( cohesion .eqv. .true. ) then
                if ( tstep .eq. 1 ) then
                    if ( deltan(j, i) .ge. -5d-1 ) then ! can be fancier
                        bond (j, i) = 1
                    end if
                    call bond_properties (j, i)
                end if
			end if

            ! verify if two particles are colliding
            if ( deltan(j,i) .gt. 0 ) then
                
                call contact_forces (j, i)
				!call bond_creation (j, i) ! to implement

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

            else
            
                call reset_contact (j, i)

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

            else

                call reset_bond (j, i)

			end if

			! compute sheltering height for particule j on particle i for air and water drag
            ! you have to check both sides of the matrix because it is not symmetric
            call sheltering(j, i)
            call sheltering(i, j)

        end do

        ! compute the total forcing from winds, currents and coriolis on particule i
        call forcing (i)
        call coriolis(i)

    end do
    !$omp end parallel do

    ! deallocate tree memory
    call tree%deallocate()

    ! reduce all the force variables
    call force_reduction

    ! sum all forces together on particule i
    do i = last_iter, first_iter, -1
        tfx_r(i) = fcx_r(i) + fbx_r(i) + fax(i) + fwx(i) + fcorx(i)
        tfy_r(i) = fcy_r(i) + fby_r(i) + fay(i) + fwy(i) + fcory(i)

        ! sum all moments on particule i together
        m_r(i) =  mc_r(i) + mb_r(i) + ma(i) + mw(i)
    end do

    ! broadcast forces to all so that they can each update their x and u
    call broadcast_total_forces

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


    !tfx(2) = tfx(2) - 1d8

    !tfx(44) = 0d0
    !tfy(44) = 0d0
    !tfx(45) = 0d0
    !tfy(45) = 0d0

    ! bottom plate on the left
    !tfy(1) = tfy(1) + 1d8
    !tfy(2) = tfy(2) + 1d8
    !tfy(3) = tfy(3) + 1d8

    ! top plate on the right
    !tfy(6) = tfy(6) - 5d7
    !tfy(4) = tfy(4) - 5d7
    !tfy(7) = tfy(7) - 5d7

    ! left plate
    !tfx(2) = tfx(2) + 1d8
    !tfx(8) = tfx(8) + 1d8
    !tfx(14) = tfx(14) + 1d8
    !tfx(20) = tfx(20) + 1d8
    !tfx(26) = tfx(26) + 1d8

    ! right plate
    !tfx(7) = tfx(7) - 5d7
    !tfx(13) = tfx(13) - 5d7
    !tfx(19) = tfx(19) - 5d7
    !tfx(25) = tfx(25) - 5d7
    !tfx(31) = tfx(31) - 5d7

end subroutine experiment_forces
