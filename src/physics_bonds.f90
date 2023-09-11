subroutine bond_properties (j, i)

	implicit none

	include "parameter.h"
	include "CB_variables.h"
	include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: i, j
    double precision :: gb

    gb = eb / 2d0 / (1 + poiss_ratio)

	rb  (j, i) = lambda_rb * min(r(i), r(j))
	hb  (j, i) = (h(i) + h(j)) / 2d0
	lb  (j, i) = lambda_lb * (r(i) + r(j))

    sb  (j, i) = 2d0 * rb (j, i) * hb (j, i)
	ib  (j, i) = 2d0 / 3d0 * hb (j, i) * rb (j, i) ** 3d0

	knb (j, i) = eb / lb (j, i)
	ktb (j, i) = 5d0 / 6d0 * gb / lb (j, i)

end subroutine bond_properties