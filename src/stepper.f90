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
    include "CB_diagnostics.h"

    integer :: i, j, k
    integer, intent(in) :: tstep
    type(KdTree) :: tree
    type(KdTreeSearch) :: search
    type(dArgDynamicArray) :: da

    ! Build the tree
    tree = KdTree(x, y)

    ! reset the forces and sheltering height
    call reset_forces
    if ( shelter .eqv. .true. ) then
        call reset_shelter
    end if
    
    ! put yourself in the referential of the ith particle
	! loop through all j particles and compute interactions

    !$omp parallel do schedule(dynamic, 1) &
    !$omp private(i,j,da) &
    !$omp reduction(+:fcx,fcy,mc,fbx,fby,mb) &
    !&&#ifdef DIAG
    !$omp reduction(+:sigxx,sigyy,sigxy,sigyx) &
    !$omp reduction(+:tac,tab,pc,pb)
    !&&#endif
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
                    if ( deltan(j, i) .ge. -1d2 ) then ! can be fancier
                        bond (j, i) = 1
                    end if
                    call bond_properties (j, i)
                end if
			end if

            ! verify if two particles are colliding
            if ( deltan(j,i) .gt. 0 ) then
                
                call contact_forces (j, i)
				!call bond_creation (j, i) ! to implement
                
                ! change coordinate system
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
                mc(j) = mc(j) - r(j) * fct(j,i) - mcc(j,i)

            else
            
                call reset_contact (j, i)

            end if

			! compute forces from bonds between particle i and j
			if ( bond (j, i) .eq. 1 ) then

				call bond_forces (j, i)
				call bond_breaking (j, i)

                if ( bond (j, i) .eq. 1 ) then

                    ! change coordinate system
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
                    mb(j) = mb(j) - r(j) * fbt(j,i) - mbb(j, i)

                end if

            else

                call reset_bond (j, i)

			end if

			! compute sheltering height for particule j on particle i for air and water drag
            ! you have to check both sides of the matrix because it is not symmetric
            if ( shelter .eqv. .true. ) then
                call sheltering(j, i)
            end if

            !-------------------------------------------------------
            !           Computation of diagnotics variables
            !-------------------------------------------------------
            
            ! if ( flag_diag_stress .eqv. .true. ) then
            ! compute the stress using cauchy stress formula (this needs to be averaged over the size of the particle)
            !-------------------------------------------------------
            !
            !    \sigma_{ij} = 1/A \sum_{c} r_j * Fcn_i
            !
            !-------------------------------------------------------

            sigxx(i) = sigxx(i) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        cosa(j,i) * r(i) * cosa(j,i)
            sigyy(i) = sigyy(i) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        sina(j,i) * r(i) * sina(j,i)
            sigxy(i) = sigxy(i) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        sina(j,i) * r(i) * cosa(j,i)
            sigyx(i) = sigyx(i) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        cosa(j,i) * r(i) * sina(j,i)

            ! Newton's third law equivalent for stress
            sigxx(j) = sigxx(j) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        cosa(j,i) * r(i) * cosa(j,i)
            sigyy(j) = sigyy(j) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        sina(j,i) * r(i) * sina(j,i)
            sigxy(j) = sigxy(j) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        sina(j,i) * r(i) * cosa(j,i)
            sigyx(j) = sigyx(j) - (                            &
                        sqrt(fcn(j,i) ** 2 + fct(j,i) ** 2) +  &
                        sqrt(fbn(j,i) ** 2 + fbt(j,i) ** 2)) * &
                        cosa(j,i) * r(i) * sina(j,i)
            ! end if
            
            ! if ( flag_diag_pressure .eqv. .true. ) then
            !     call diag_mean_pressure (j, i)
            ! end if

        end do

        ! compute the total forcing from winds, currents and coriolis on particule i
        call forcing (i)
        call coriolis(i)

         ! verify the bondary conditions for each particle
        call verify_bc (i)

    end do
    !$omp end parallel do

    ! deallocate tree memory
    call tree%deallocate()

    ! reduce all the force variables
    call force_reduction

    ! sum all forces together on particule i
    do i = last_iter, first_iter, -1
        tfx_r(i) = fcx_r(i) + fbx_r(i) + fax(i) + fwx(i) + fcorx(i) &
                    + fx_bc(i)
        tfy_r(i) = fcy_r(i) + fby_r(i) + fay(i) + fwy(i) + fcory(i) &
                    + fy_bc(i)

        ! sum all moments on particule i together
        m_r(i) =  mc_r(i) + mb_r(i) + ma(i) + mw(i) + m_bc(i)

        ! same for stresses
        tsigxx_r(i) = sigxx_r(i) + sigxx_bc(i)
        tsigyy_r(i) = sigyy_r(i) + sigyy_bc(i)
        tsigxy_r(i) = sigxy_r(i) + sigxy_bc(i)
        tsigyx_r(i) = sigyx_r(i) + sigyx_bc(i)

        ! ! same for pressure
        ! tp_r(i) = pc_r(i) + pc_r(i) + p_bc(i)
    end do

    ! broadcast forces to all so that the nodes can each update their x and u
    call broadcast_total_forces

    ! forces on right side plate
    call normal_forces("right", tstep)

    ! forces on left side plate
    !call normal_forces("left", tstep)

    ! integration in time
    call velocity
    call position

end subroutine stepper


subroutine normal_forces (side, tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"

    character (*), intent(in) :: side
    integer, intent(in) :: tstep

    double precision :: ftmp, tau
    integer :: i
    
    ! ~12 seconds is required to resolve the elastic wave travelling 30km
    tau = 50 / dt

    ! forces on the right of the assembly
    ! forces on the right of the assembly
    if ( side == "right" ) then
        ftmp = 0d0!maxval(tfx(n-29:n))
!        if (tstep*dt < 75d-2) then
!            if (tstep*dt > 25d-2) then
                !do i = n, n - 9, -1!n, n - 29 , -1
                    tfx(10) = tfx(10) + ftmp + pfn * tanh( tstep / tau )
                !end do
!            end if       
        !end if
     end if    

    ! no forces on the left of assembly (fixed)
    if ( side == "left" ) then
        do i = n - 30, n - 59 , -1
            tfy(i) = 0d0
            tfx(i) = 0d0
        end do
    end if


end subroutine normal_forces