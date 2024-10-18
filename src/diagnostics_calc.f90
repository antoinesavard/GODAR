subroutine diag_stress (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: j, i
    
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

end subroutine diag_stress


subroutine diag_mean_pressure (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_diagnostics.h"

    integer, intent(in) :: j, i

    ! compute the average pressure inside particle i
    !-------------------------------------------------------
    !
    !    P_i = \sum_{c} Fcn_{ij} * a_{ij} / \sum_{c} a_{ij}
    !
    !-------------------------------------------------------
    ac     = deltat(j, i) * min(h(i), h(j))
    tac(i) = tac(i) + ac

    tab(i) = tab(i) + sb(j, i)
    
    pc(i)  = pc(i) - fcn(j, i) * ac
    pb(i)  = pb(i) - fbn(j, i) * sb(j, i)

    tp(i)  = tp(i) + pc(i) / tac(i)        &
                   + pb(i) / tab(i)       

    ! symmetric part
    tac(j) = tac(j) + ac

    tab(j) = tab(j) + sb(j, i)
    
    pc(j)  = pc(j) - fcn(j, i) * ac
    pb(j)  = pb(j) - fbn(j, i) * sb(j, i)

    tp(j)  = tp(j) + pc(i) / tac(i)        &
                   + pb(i) / tab(i)


end subroutine diag_mean_pressure
