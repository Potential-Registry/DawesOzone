module splin
    ! use path_resolution_mod
    implicit none

    ! Global variables within this module
    integer :: splin_nr, splin_nlt   ! base points for interpolation
!    integer, parameter :: splin_nlt = 45  ! nombre de  coefficients ds tabCn-roo.res
    parameter(splin_nr=16, splin_nlt=45)

CONTAINS

    subroutine inpdat(coef, ll, qq2, rro)
        !Initialization of the input data
        implicit none
        !integer, parameter :: nr=16      ! base points for interpolation
        integer i, nl, n
        integer :: file_unit
        !double precision, dimension (nr) :: rro(nr), qq2(nr)
        double precision yy, zz, rmin, rmax, rstep


        double precision :: coef(splin_nlt, splin_nr)
        integer :: ll(splin_nlt, 6)
        double precision xmin, xmax
        double precision :: rro(splin_nr), qq2(splin_nr)

        ! compiled in data for speed and path-independence

        ll = reshape((/&
              -1,-1,1,0,0,1,&
              -1,-1,1,0,0,2,&
              -1,0,1,-1,0,1,&
              -1,0,1,-1,0,2,&
              -1,1,1,-2,0,1,&
              -1,1,1,-2,0,2,&
              0,-1,1,1,0,1,&
              0,-1,1,1,0,2,&
              0,0,1,0,0,1,&
              0,0,1,0,0,2,&
              0,1,1,-1,0,1,&
              0,1,1,-1,0,2,&
              1,-1,1,2,0,1,&
              1,-1,1,2,0,2,&
              1,0,1,1,0,1,&
              1,0,1,1,0,2,&
              1,1,1,0,0,1,&
              1,1,1,0,0,2,&
              -1,-1,1,0,0,0,&
              -1,-1,1,0,0,1,&
              -1,-1,1,0,0,2,&
              -1,0,1,-1,0,0,&
              -1,0,1,-1,0,1,&
              -1,0,1,-1,0,2,&
              -1,1,1,-2,0,0,&
              -1,1,1,-2,0,1,&
              -1,1,1,-2,0,2,&
              0,-1,1,1,0,0,&
              0,-1,1,1,0,1,&
              0,-1,1,1,0,2,&
              0,0,1,0,0,0,&
              0,0,1,0,0,1,&
              0,0,1,0,0,2,&
              0,1,1,-1,0,0,&
              0,1,1,-1,0,1,&
              0,1,1,-1,0,2,&
              1,-1,1,2,0,0,&
              1,-1,1,2,0,1,&
              1,-1,1,2,0,2,&
              1,0,1,1,0,0,&
              1,0,1,1,0,1,&
              1,0,1,1,0,2,&
              1,1,1,0,0,0,&
              1,1,1,0,0,1,&
              1,1,1,0,0,2/),&
          shape(ll))
        coef= reshape((/&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -2.88086,-2.43796,-1.99981,-1.56133,-1.12175,-0.682985,-0.248346,0.178082,0.592202,0.990201,1.36869,1.72489,2.05632,2.36071,2.6357,2.87846,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              3.32653,2.81511,2.30918,1.80287,1.29529,0.788644,0.286765,-0.205632,-0.683815,-1.14339,-1.58043,-1.99173,-2.37443,-2.72591,-3.04345,-3.32376,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -1.1761,-0.995292,-0.816419,-0.637409,-0.457953,-0.278828,-0.101387,0.0727018,0.241765,0.404248,0.558766,0.704182,0.839488,0.963754,1.07602,1.17513,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -3.32653,-2.81511,-2.30918,-1.80287,-1.29529,-0.788644,-0.286765,0.205632,0.683815,1.14339,1.58043,1.99173,2.37443,2.72591,3.04345,3.32376,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              5.76171,4.87592,3.99962,3.12265,2.2435,1.36597,0.496692,-0.356164,-1.1844,-1.9804,-2.73739,-3.44977,-4.11264,-4.72141,-5.27141,-5.75693,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -3.32653,-2.81511,-2.30918,-1.80287,-1.29529,-0.788644,-0.286765,0.205632,0.683815,1.14339,1.58043,1.99173,2.37443,2.72591,3.04345,3.32376,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -1.1761,-0.995292,-0.816419,-0.637409,-0.457953,-0.278828,-0.101387,0.0727018,0.241765,0.404248,0.558766,0.704182,0.839488,0.963754,1.07602,1.17513,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              3.32653,2.81511,2.30918,1.80287,1.29529,0.788644,0.286765,-0.205632,-0.683815,-1.14339,-1.58043,-1.99173,-2.37443,-2.72591,-3.04345,-3.32376,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -2.88086,-2.43796,-1.99981,-1.56133,-1.12175,-0.682985,-0.248346,0.178082,0.592202,0.990201,1.36869,1.72489,2.05632,2.36071,2.6357,2.87846,&
              -26.8054,-27.2254,-27.8935,-28.7175,-29.6734,-30.4767,-31.4433,-32.3389,-33.4974,-34.1125,-34.9134,-33.8497,-30.7686,-36.2997,-36.7577,-36.9238,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -3.0483,-3.26079,-3.54863,-3.88964,-4.2206,-4.65401,-5.07147,-5.47235,-5.7317,-6.33899,-6.74596,-5.35107,-2.06491,-7.37353,-7.44939,-7.56679,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -0.18057,-0.19306,-0.2116,-0.234275,-0.259626,-0.28703,-0.315466,-0.344297,-0.370949,-0.4014,-0.42969,-0.441283,-0.429574,-0.495326,-0.509486,-0.519575,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.127682,0.136514,0.149624,0.165658,0.183583,0.20296,0.223068,0.243455,0.262301,0.283832,0.303837,0.312034,0.303755,0.350248,0.360261,0.367395,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.18057,0.19306,0.2116,0.234275,0.259626,0.28703,0.315466,0.344297,0.370949,0.4014,0.42969,0.441283,0.429574,0.495326,0.509486,0.519575,&
              -27.9101,-28.3449,-29.0429,-29.9062,-30.9069,-31.7577,-32.7737,-33.7184,-34.9288,-35.586,-36.4313,-35.3954,-32.3193,-37.9206,-38.4033,-38.5866,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -3.51743,-3.76237,-4.09838,-4.49831,-4.89512,-5.39973,-5.89108,-6.36686,-6.69545,-7.38186,-7.86233,-6.49756,-3.18098,-8.66042,-8.77308,-8.91668,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.18057,0.19306,0.2116,0.234275,0.259626,0.28703,0.315466,0.344297,0.370949,0.4014,0.42969,0.441283,0.429574,0.495326,0.509486,0.519575,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.127682,0.136514,0.149624,0.165658,0.183583,0.20296,0.223068,0.243455,0.262301,0.283832,0.303837,0.312034,0.303755,0.350248,0.360261,0.367395,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -0.18057,-0.19306,-0.2116,-0.234275,-0.259626,-0.28703,-0.315466,-0.344297,-0.370949,-0.4014,-0.42969,-0.441283,-0.429574,-0.495326,-0.509486,-0.519575,&
              -26.8054,-27.2254,-27.8935,-28.7175,-29.6734,-30.4767,-31.4433,-32.3389,-33.4974,-34.1125,-34.9134,-33.8497,-30.7686,-36.2997,-36.7577,-36.9238,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
              -3.0483,-3.26079,-3.54863,-3.88964,-4.2206,-4.65401,-5.07147,-5.47235,-5.7317,-6.33899,-6.74596,-5.35107,-2.06491,-7.37353,-7.44939,-7.56679/),&
          shape(coef))

!        read(file_unit, *)
!        do nl = 1, nlt
!            read(file_unit, *)(ll(nl, i), i = 1, 6), (coef(nl, n), n = 1, nr)
!        enddo
!        close(file_unit)

        rmin = 1.8
        rmax = 3.3
        rstep = (rmax - rmin) / splin_nr
        do i=1, splin_nr
          rro(i) = rmin + rstep * (i-1)
        end do
        qq2 = (/-1.01083,-0.855424,-0.701688,-0.547834,-0.393597,-0.239644,-0.087139,0.062485,0.20779,0.347439,0.480243,0.605223,0.721515,0.828318,0.924808,1.00999/)

    end subroutine inpdat
    !
    subroutine initspline(np, xi, yi, b, c, d, coef)
        ! Initialization of the spline
        implicit none

        double precision:: xi(splin_nr), yi(splin_nr), b(splin_nr), c(splin_nr), d(splin_nr)
        double precision :: xmin, xmax, step
        double precision :: coef(splin_nlt, splin_nr)
        integer :: i, np

        xmin = 1.8
        xmax = 3.3

        step = (xmax - xmin) / (splin_nr - 1)
        do i = 1, splin_nr
            xi(i) = xmin + step * float(i - 1)
            yi(i) = coef(np, i)
            !  write (*,200) xi(i), yi(i)
        end do
        !  step 2: call spline to calculate spline coeficients
        call spline_d (xi, yi, b, c, d, splin_nr)
!        200 format (3f12.5)
    end subroutine initspline
    !
    subroutine spline_d (x, y, b, c, d, n)
        !======================================================================
        !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
        !  for cubic spline interpolation
        !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
        !  for  x(i) <= x <= x(i+1)
        !  Alex G: January 2010
        !----------------------------------------------------------------------
        !  input..
        !  x = the arrays of data abscissas (in strictly increasing order)
        !  y = the arrays of data ordinates
        !  n = size of the arrays xi() and yi() (n>=2)
        !  output..
        !  b, c, d  = arrays of spline coefficients
        !  comments ...
        !  spline.f90 program is based on fortran version of program spline.f
        !  the accompanying function fspline can be used for interpolation
        !======================================================================
        implicit none
        integer n
        double precision x(n), y(n), b(n), c(n), d(n)
        integer i, j, gap
        double precision h

        gap = n - 1
        ! check input
        if (n < 2) return
        if (n < 3) then
            b(1) = (y(2) - y(1)) / (x(2) - x(1))   ! linear interpolation
            c(1) = 0.
            d(1) = 0.
            b(2) = b(1)
            c(2) = 0.
            d(2) = 0.
            return
        end if
        !
        ! step 1: preparation
        !
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1)) / d(1)
        do i = 2, gap
            d(i) = x(i + 1) - x(i)
            b(i) = 2.0 * (d(i - 1) + d(i))
            c(i + 1) = (y(i + 1) - y(i)) / d(i)
            c(i) = c(i + 1) - c(i)
        end do
        !
        ! step 2: end conditions
        !
        b(1) = -d(1)
        b(n) = -d(n - 1)
        c(1) = 0.0
        c(n) = 0.0
        if(n /= 3) then
            c(1) = c(3) / (x(4) - x(2)) - c(2) / (x(3) - x(1))
            c(n) = c(n - 1) / (x(n) - x(n - 2)) - c(n - 2) / (x(n - 1) - x(n - 3))
            c(1) = c(1) * d(1)**2 / (x(4) - x(1))
            c(n) = -c(n) * d(n - 1)**2 / (x(n) - x(n - 3))
        end if
        !
        ! step 3: forward elimination
        !
        do i = 2, n
            h = d(i - 1) / b(i - 1)
            b(i) = b(i) - h * d(i - 1)
            c(i) = c(i) - h * c(i - 1)
        end do
        !
        ! step 4: back substitution
        !
        c(n) = c(n) / b(n)
        do j = 1, gap
            i = n - j
            c(i) = (c(i) - d(i) * c(i + 1)) / b(i)
        end do
        !
        ! step 5: compute spline coefficients
        !
        b(n) = (y(n) - y(gap)) / d(gap) + d(gap) * (c(gap) + 2.0 * c(n))
        do i = 1, gap
            b(i) = (y(i + 1) - y(i)) / d(i) - d(i) * (c(i + 1) + 2.0 * c(i))
            d(i) = (c(i + 1) - c(i)) / d(i)
            c(i) = 3. * c(i)
        end do
        c(n) = 3.0 * c(n)
        d(n) = d(n - 1)
    end subroutine spline_d

    function ispline(u, x, y, b, c, d, n)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        implicit none
        double precision ispline
        integer n
        double precision  u, x(n), y(n), b(n), c(n), d(n)
        integer i, j, k
        double precision dx

        ! if u is ouside the x() interval take a boundary value (left or right)
        if(u <= x(1)) then
            ispline = y(1)
            return
        end if
        if(u >= x(n)) then
            ispline = y(n)
            return
        end if

        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n + 1
        do while (j > i + 1)
            k = (i + j) / 2
            if(u < x(k)) then
                j = k
            else
                i = k
            end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = u - x(i)
        ispline = y(i) + dx * (b(i) + dx * (c(i) + dx * d(i)))
    end function ispline
    !
    subroutine initspline2(xi, yi, b, c, d, rro)
        ! Initialization of the spline
        implicit none
        !integer, parameter :: nr=16      ! base points for interpolation
        double precision:: xi(splin_nr), yi(splin_nr), b(splin_nr), c(splin_nr), d(splin_nr)
        !double precision, dimension (nr) :: rro(nr),  qq2(nr)
        double precision step, xmin, xmax
        double precision :: rro(splin_nr), qq2(splin_nr)
        integer i

        xmin = 1.8
        xmax = 3.3

        step = (xmax - xmin) / (splin_nr - 1)
        do i = 1, splin_nr
            xi(i) = xmin + step * float(i - 1)
            if(dabs(xi(i) - rro(i))>=1.e-01)then
                print *, '# Erreur de distance O-O', xi(i), rro(i)
                stop
            endif
            yi(i) = qq2(i)
            !  write (*,200) xi(i), yi(i)
        end do
        !  step 2: call spline to calculate spline coeficients
        call spline_d (xi, yi, b, c, d, splin_nr)
!        200 format (3f12.5)
    end subroutine initspline2
end module splin
