subroutine diag_stress (j, i, p)

    use parameters
    use variables
    use const
    use bonds
    use diagnostics
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: j, i
    type(pair_t), intent(in) :: p
    
    ! compute the stress using cauchy stress formula (this needs to be averaged over the size of the particle)
    !-------------------------------------------------------
    !
    !    \sigma_{ij} = 1/A \sum_{c} r_j * Fcn_i
    !
    !-------------------------------------------------------

    sigxx(i) = sigxx(i) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%cosa * r(i) * p%cosa
    sigyy(i) = sigyy(i) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%sina * r(i) * p%sina
    sigxy(i) = sigxy(i) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%sina * r(i) * p%cosa
    sigyx(i) = sigyx(i) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%cosa * r(i) * p%sina

    ! Newton's third law equivalent for stress
    sigxx(j) = sigxx(j) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%cosa * r(i) * p%cosa
    sigyy(j) = sigyy(j) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%sina * r(i) * p%sina
    sigxy(j) = sigxy(j) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%sina * r(i) * p%cosa
    sigyx(j) = sigyx(j) - (                            &
                sqrt(p%fcn ** 2 + p%fct ** 2) +        &
                sqrt(p%fbn ** 2 + p%fbt ** 2)) *       &
                p%cosa * r(i) * p%sina

end subroutine diag_stress


subroutine diag_mean_pressure (j, i, p)

    use parameters
    use variables
    use const
    use bonds
    use diagnostics
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: j, i
    type(pair_t), intent(in) :: p

    double precision :: ac_ij

    ! compute the average pressure inside particle i
    !-------------------------------------------------------
    !
    !    P_i = \sum_{c} Fcn_{ij} * a_{ij} / \sum_{c} a_{ij}
    !
    !-------------------------------------------------------
    ! local area
    ac_ij = p%delt_ridge * min(h(i), h(j))
    ! total contact area
    tac(i)  = tac(i) + ac_ij
    ! total bond area
    tab(i)  = tab(i) + sb(j, i)
    
    ! pressure from contacts and bonds
    pc(i)   = pc(i) - p%fcn * ac_ij
    pb(i)   = pb(i) - p%fbn * sb(j, i)

    ! symmetric part
    tac(j) = tac(j) + ac_ij

    tab(j) = tab(j) + sb(j, i)
    
    pc(j)  = pc(j) - p%fcn * ac_ij
    pb(j)  = pb(j) - p%fbn * sb(j, i)
    ! end if


end subroutine diag_mean_pressure


subroutine dilation (j, i)

    use parameters
    use variables
    use const
    use bonds
    use diagnostics

    implicit none


    integer, intent(in) :: i, j

    

end subroutine dilation
