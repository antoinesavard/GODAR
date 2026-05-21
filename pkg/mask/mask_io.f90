module mask_io

    implicit none

    integer :: nx_mask = 0, ny_mask = 0
    double precision :: dx_mask, x_origin, y_origin, Lx, Ly
    character(len=8) :: mask_proj

    real(kind=4), allocatable :: sdf(:,:)

    double precision, parameter :: sdf_outside = huge(1.0d0)

contains

    subroutine load_sdf(filename)

        implicit none

        character(len=*), intent(in) :: filename

        integer :: u, ix, iy
        double precision :: ymax_in

        open(newunit=u, file=filename, status='old', action='read')

        read(u, *) nx_mask, ny_mask, dx_mask, x_origin, ymax_in, mask_proj

        Lx       = nx_mask * dx_mask
        Ly       = ny_mask * dx_mask
        y_origin = ymax_in - Ly

        allocate(sdf(nx_mask, ny_mask))

        ! file rows are top-to-bottom, flip so iy=1 is at the bottom
        do iy = ny_mask, 1, -1
            read(u, *) (sdf(ix, iy), ix = 1, nx_mask)
        end do

        close(u)

    end subroutine load_sdf


    pure function sdf_at(x, y) result(d)

        ! Bilinear interp of sdf, treating values as samples at cell
        ! centers. Outside the grid extent: returns sdf_outside.

        double precision, intent(in) :: x, y
        double precision :: d

        double precision :: fx, fy, tx, ty
        integer :: ix0, iy0, ix1, iy1

        if (x < x_origin .or. x > x_origin + nx_mask * dx_mask .or. &
            y < y_origin .or. y > y_origin + ny_mask * dx_mask) then
            d = sdf_outside
            return
        end if

        fx = (x - x_origin) / dx_mask - 0.5d0
        fy = (y - y_origin) / dx_mask - 0.5d0

        ix0 = max(1, min(nx_mask - 1, floor(fx) + 1))
        iy0 = max(1, min(ny_mask - 1, floor(fy) + 1))
        ix1 = ix0 + 1
        iy1 = iy0 + 1

        tx = max(0.0d0, min(1.0d0, fx - dble(ix0 - 1)))
        ty = max(0.0d0, min(1.0d0, fy - dble(iy0 - 1)))

        d = (1.0d0 - tx) * (1.0d0 - ty) * dble(sdf(ix0, iy0)) &
          +          tx  * (1.0d0 - ty) * dble(sdf(ix1, iy0)) &
          + (1.0d0 - tx) *          ty  * dble(sdf(ix0, iy1)) &
          +          tx  *          ty  * dble(sdf(ix1, iy1))

    end function sdf_at


    pure subroutine sdf_grad(x, y, gx, gy)

        double precision, intent(in)  :: x, y
        double precision, intent(out) :: gx, gy

        double precision :: h

        h  = dx_mask
        gx = (sdf_at(x + h, y) - sdf_at(x - h, y)) / (2.0d0 * h)
        gy = (sdf_at(x, y + h) - sdf_at(x, y - h)) / (2.0d0 * h)

    end subroutine sdf_grad


    pure function lat_at(x, y) result(phi)

        ! North Pole at (0,0) for both EASE2-N and PS-N (NSIDC).
        ! Returns latitude in radians on a spherical Earth (R = 6371 km).

        double precision, intent(in) :: x, y
        double precision :: phi

        double precision :: r, ratio, m_pole
        double precision, parameter :: R_earth   = 6371000.0d0
        double precision, parameter :: pi_half   = 1.5707963267948966d0
        double precision, parameter :: ps_truelat = 70.0d0 * 3.141592653589793d0 / 180.0d0

        r = sqrt(x * x + y * y)

        select case (trim(mask_proj))
        case ("EASE2-N")
            ratio = min(1.0d0, r / (2.0d0 * R_earth))
            phi   = pi_half - 2.0d0 * asin(ratio)
        case ("PS-N")
            m_pole = (1.0d0 + sin(ps_truelat)) / 2.0d0
            phi    = pi_half - 2.0d0 * atan2(r, 2.0d0 * R_earth * m_pole)
        case default
            ! fall back to namelist lat
            phi = -huge(1.0d0)
        end select

    end function lat_at

end module mask_io
