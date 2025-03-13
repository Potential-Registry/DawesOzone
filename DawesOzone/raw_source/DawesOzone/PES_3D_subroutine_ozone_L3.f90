!MODULE dynamic_parameters
!    implicit none
!    !save
!    public
!    real(kind = 8), allocatable :: b2(:, :), b2_lower(:, :), b2_minimal(:, :), b2_seed(:, :), d(:), d_seed(:), Jac(:), Jac2(:), coords(:, :), coords_seed(:, :)
!    real(kind = 8), allocatable :: cart(:), bdist(:), stored_weights(:), stored_grad(:, :)
!    real(kind = 8), allocatable :: pot(:), pot_seed(:), grad(:, :), grad_seed(:, :), dcart(:)
!    integer :: basis_1, basis_2, basis_3, order_1, order_2, order_3, count3, zz, zz_low, zz4, support, count7, myid, natom, ab_flag, ab_flag2, lab, permfac
!    integer :: maxpoints, nbdist, order_1_min, order_2_min, order_3_min, count_seed, low_grid, subzero
!    real(kind = 8) :: epss, acc, rmax(3), rmin(3), poten, hartokcl, ugrad, Max_E, Max_E_seed, E_range, dist_tol(3), E_limit, alpha(3), Glob_min, current_geom(3), current_geom_grad(3), mass(3)
!    character(len = 3) :: symb(3)
!    integer, allocatable :: stored_weights_ind(:)
!END MODULE dynamic_parameters

! the actual function that does the potential calls
! for thread-safety/cleanliness I eliminated the global variables
! and instead set up two subroutines that handle the data that IMLS reads in
subroutine IMLS(jac3, V, SO_flag)
    implicit none

    integer, intent(in) :: SO_flag
    double precision, dimension(3), intent(in) :: jac3
    double precision, intent(out) :: V

    integer :: i, initflag, data_unit
    real(kind=8) :: tocm

    character(len = 27) :: path

    real(kind = 8) :: epss, pii, hartokcl, alpha(3)
    integer :: order_1_min, order_2_min, order_3_min

    integer :: count3, order_1, order_2, order_3, basis_1, basis_2, basis_3, maxpoints, lab, zz, low_grid, count_seed
    real(kind = 8) :: rmax(3), rmin(3), mass(3), Max_E

!    integer :: count32, order_12, order_22, order_32, basis_12, basis_22, basis_32, maxpoints2, lab2, zz2, low_grid2, count_seed2
!    real(kind = 8) :: rmax2(3), rmin2(3), mass2(3), Max_E2

!    logical data_open

!    path = 'data/PES_2663_QZ_F12_20.dat'
!    data_unit = 652
!    inquire(unit=data_unit, opened=data_open)
!    if (data_open) then
!        data_unit = data_unit + irand()
!    end if
!    inquire(unit=data_unit, opened=data_open)
!    if (data_open) then
!        print*, "........?", V
!        error stop "wat data lready open???"
!    end if

    count3 = 5326
    order_1 = 12
    order_2 = 9
    order_3 = 7
    maxpoints = 7000
    lab = 2
    zz = 4
    low_grid = 0
    count_seed = 0
    mass = (/15.994914622100000, 15.994914622100000, 15.994914622100000/)
    rmax = (/ 1.7000000000000000, 10.000000000000000, 0.50000000000000000/)
    rmin = (/ 1.0000000000000000, 1.0000000000000000, -1.0000000000000000/)

!    OPEN(UNIT = data_unit, FILE = path, FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL')
!    read(data_unit) count32
!    read(data_unit) order_12
!    read(data_unit) order_22
!    read(data_unit) order_32
!    read(data_unit) maxpoints2
!    read(data_unit) mass2
!    read(data_unit) rmax2
!    read(data_unit) rmin2
!    read(data_unit) Max_E2
!    read(data_unit) lab2
!    read(data_unit) zz2
!    read(data_unit) low_grid2
!    read(data_unit) count_seed2

!    print *, count3, order_1, order_2, order_3, maxpoints
!    print *, lab, zz, low_grid, count_seed
!    print *, mass
!    print *, rmax
!    print *, rmin

    order_1_min = 3  !!!!!basis used to fit low level grid.
    order_2_min = 3
    order_3_min = 3

    call basis_size(3, order_1, order_2, order_3, lab, basis_1)
    call basis_size(3, order_1 - 1, order_2 - 1, order_3 - 1, lab, basis_2)
    call basis_size(3, order_1_min, order_2_min, order_3_min, lab, basis_3)

!    print *, basis_1, count3, basis_2, basis_3
!
!    error stop "wat"

    ! IMLS_driver isn't really a different subroutine, but rather
    ! the second half of this one, but with direct knowledge of array dimensions and stuff
    ! it then calls IMLS_core which _is_ conceptually different, since it actively gets the
    ! potential for the loaded data
    call IMLS_driver(&
            jac3, SO_flag, V, data_unit, &
            count3, order_1, order_2, order_3, &
            order_1_min, order_2_min, order_3_min, &
            basis_1, basis_2, basis_3,&
            maxpoints, count_seed, mass, rmax, rmin, Max_E, lab, zz, low_grid)

!    print*, ">>>>", V

    contains

        subroutine IMLS_driver(&
            jac3, SO_flag, V, data_unit, &
            count3, order_1, order_2, order_3, &
            order_1_min, order_2_min, order_3_min, &
            basis_1, basis_2, basis_3,&
            maxpoints, count_seed, mass, rmax, rmin, Max_E, lab, zz, low_grid)
        implicit none

        integer, intent(in) :: SO_flag, data_unit, count3, basis_1, basis_2, basis_3, &
                order_1, order_2, order_3, order_1_min, order_2_min, order_3_min, &
                maxpoints, count_seed, lab, zz, low_grid
        real(kind=8), intent(in) :: jac3(:), mass(:), rmax(:), rmin(:), Max_E
        real(kind=8), intent(out) :: V

        real(kind = 8) :: b2(basis_1, count3), b2_lower(basis_2, count3), b2_minimal(basis_3, count3), d(count3), coords(count3, 3)
        real(kind = 8) :: b2_seed(basis_3, count_seed), d_seed(count_seed), coords_seed(count_seed, 3), Max_E_seed

        logical data_open

        b2_lower = 0d0
        b2_minimal = 0d0
        b2_seed = 0d0
        d = 0d0
        d_seed = 0d0
        coords = 0d0
        coords_seed = 0d0

!        do i=1,count3
!          read(data_unit) b2(:,i)
!        enddo
!        do i=1,count3
!          read(data_unit) b2_lower(:,i)
!        enddo
!        do i=1,count3
!          read(data_unit) b2_minimal(:,i)
!        enddo
!        do i=1,count3
!          read(data_unit) d(i)
!        enddo
!        do i=1,count3
!          read(data_unit) coords(i,:)
!     !      write(223,'(4f20.15)') coords(i,:)
!        enddo
     !stop
!        if(low_grid>0)then
!          read(data_unit) Max_E_seed
!          do i=1,count_seed
!            read(data_unit) b2_seed(:,i)
!          enddo
!          do i=1,count_seed
!            read(data_unit) d_seed(i)
!          enddo
!          do i=1,count_seed
!            read(data_unit) coords_seed(i,:)
!          enddo
!        endif
!        close(data_unit)
!
!        print*, ">oh?", b2(1, :10)
!        print*, ">oh.?", b2_lower(10, :10)
!        print*, ">oh...?", b2_minimal(10, :10)
!        print*, ">d...?", d(count3-5:count3)
!        print*, ">crds...2?", coords(count3-5:count3, :)

        ! fills in fit data
        call fill_b2(b2, basis_1, count3)
        call fill_b2_lower(b2_lower, basis_2, count3)
        call fill_b2_minimal(b2_minimal, basis_3, count3)
        call fill_d(d, count3)
        call fill_coords(coords, count3)

!        print*, ">oh2?", b2(1, :10)
!        print*, ">oh.2?", b2_lower(10, :10)
!        print*, ">oh...2?", b2_minimal(10, :10)
!        print*, ">d...?", d(count3-5:count3)
!        print*, ">crds...2?", coords(count3-5:count3, :)

        call IMLS_core(&
            jac3, SO_flag, V,&
            count3, order_1, order_2, order_3, &
            order_1_min, order_2_min, order_3_min, &
            basis_1, basis_2, basis_3,&
            maxpoints, count_seed, mass, rmax, rmin, Max_E, lab, zz, low_grid,&
            b2, b2_lower, b2_minimal, d, coords,&
            b2_seed, d_seed, coords_seed, Max_E_seed&
            )

    end subroutine IMLS_driver

end subroutine IMLS

! take the data from the source file now that we have dimension info and
! feed it into the core
subroutine IMLS_core(&
            input_coords, SO_flag, out_pot,&
            count3, order_1, order_2, order_3, &
            order_1_min, order_2_min, order_3_min, &
            basis_1, basis_2, basis_3,&
            maxpoints, count_seed, mass, rmax, rmin, Max_E, lab, zz, low_grid,&
            b2, b2_lower, b2_minimal, d, coords,&
            b2_seed, d_seed, coords_seed, Max_E_seed&
            )
    use splin
    use poteval
    implicit none

    real(kind=8), intent(out) :: out_pot

    integer, intent(in) :: SO_flag, count3, basis_1, basis_2, basis_3, &
            order_1, order_2, order_3, order_1_min, order_2_min, order_3_min, &
            maxpoints, count_seed, lab, zz, low_grid
    real(kind=8), intent(in) :: input_coords(3), rmax(3), rmin(3), mass(3), Max_E

    real(kind = 8), intent(in) :: b2(basis_1, count3), b2_lower(basis_2, count3), b2_minimal(basis_3, count3), d(count3), coords(count3, 3)
    real(kind = 8), intent(in) :: b2_seed(basis_3, count_seed), d_seed(count_seed), coords_seed(count_seed, 3), Max_E_seed

    ! new parameters
    double precision :: xi(3), temp3, h2wn, cartt(9), dist(3), dist_temp(3), vec1(2), vec2(2), bohr, pii
    double precision :: theta(3), SO_corr, SO_corr2, cartt2(9), somme, vx, vxso, SS, tocm, rooi

    integer :: zz_low, zz4, support, count7, myid, natom, ab_flag, ab_flag2, permfac
    integer :: nbdist, subzero

    real(kind=8) :: epss, jac(3), jac2(3), alpha(3), current_geom(3), hartokcl, dist_tol(3)
    integer :: stored_weights_ind(count3)


    double precision :: theta(3), SO_corr, SO_corr2, cartt2(9), somme, vx, vxso, SS, tocm, rooi

    double precision  ::  coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
    integer :: ll(splin_nlt, 6)

    character(len = 3) :: symb(3)

    integer :: i

!    common/drob/q20O2, rooi

    call dataread
    call inpdat(coef, ll, qq2, rro, splin_nr, splin_nlt)

    current_geom = 0d0
    rooi = 0.d0
    bohr = 0.529177249d0
    ! bohr = ang2bohr
    hartokcl = 627.5095d0
    h2wn = 219474.63d0
    ! h2wn = autown
    pii = acos(-1d0)
    alpha(:) = -1d0

    epss = 1d-14
    zz_low = 4
    zz4 = 20
    natom = 3
    dist_tol = 0.7d0

    xi = input_coords

    xi(1) = xi(1) * bohr
    xi(2) = xi(2) * bohr
    xi(3) = pii * xi(3) / 180d0

    cartt = 0d0
    cartt(1) = xi(1)
    cartt(7) = xi(2) * dcos(xi(3))
    cartt(8) = xi(2) * dsin(xi(3))
    dist_temp(1) = xi(1)
    dist_temp(2) = xi(2)
    dist_temp(3) = dsqrt((cartt(1) - cartt(7))**2 + cartt(8)**2)
    theta(1) = dcos(xi(3))
    !theta(2)=((cartt(7)-cartt(4))*(-cartt(1)))/(dist_temp(1)*dist_temp(3))
    theta(2) = (dist_temp(1)**2 + dist_temp(3)**2 - dist_temp(2)**2) / (2d0 * dist_temp(1) * dist_temp(3))
    theta(3) = dcos(pii - dacos(theta(1)) - dacos(theta(2)))

    if(dist_temp(1)<=dist_temp(2))then
        if(dist_temp(1)<=dist_temp(3))then
            if(dist_temp(2)<=dist_temp(3))then
                xi(1) = dist_temp(1)
                xi(2) = dist_temp(2)
                xi(3) = theta(1)
            else
                xi(1) = dist_temp(1)
                xi(2) = dist_temp(3)
                xi(3) = theta(2)
            endif
        else
            xi(1) = dist_temp(3)
            xi(2) = dist_temp(1)
            xi(3) = theta(2)
        endif
    else
        if(dist_temp(2)<=dist_temp(3))then
            if(dist_temp(1)<=dist_temp(3))then
                xi(1) = dist_temp(2)
                xi(2) = dist_temp(1)
                xi(3) = theta(1)
            else
                xi(1) = dist_temp(2)
                xi(2) = dist_temp(3)
                xi(3) = theta(3)
            endif
        else
            xi(1) = dist_temp(3)
            xi(2) = dist_temp(2)
            xi(3) = theta(3)
        endif
    endif

    do i = 1, 1
        if(xi(i)<rmin(i)) xi(i) = rmin(i)
        if(xi(i)>rmax(i)) xi(i) = rmax(i)
    enddo

!    print*, "huh", cartt
!    print*, "huh", xi, rmin, rmax

    if(xi(2)<rmax(2)) then
        temp3 = func_actual_min(xi)

!        print*, "temp 3", temp3

        if(temp3>Max_E)then
            out_pot = ((Max_E + 141338.60419d0) / hartokcl) * h2wn
            !   write(*,*) 'hit ceiling on low grid'
            !   v=Max_E
            return
        endif
        !!!!!!!!!!!!!!!!
        temp3 = func_actual(xi)
        if(temp3>Max_E)then
            temp3 = Max_E
            !   write(*,*) 'hit ceiling'
        endif

        if(SO_flag>0)then
            out_pot = ((temp3 + 141338.60617d0) / hartokcl) * h2wn
        else
            out_pot = ((temp3 + 141338.60419d0) / hartokcl) * h2wn
        endif
    else
        out_pot = 9270.6d0
    endif


    call INT_Cart(xi, cartt, mass, 2)
    cartt2 = cartt
    cartt2(1:3) = cartt(7:9)
    cartt2(7:9) = cartt(1:3)

    call Cart_INT(xi, cartt2, mass, 1)
    xi(3) = abs(xi(3))
    xi(3) = dacos(xi(3))
    xi(3) = 180d0 * xi(3) / pii
    !write(*,*) xi
    xi(2) = xi(2) / bohr

    SO_corr = 0d0
    SO_corr2 = 0d0
    if(SO_flag>0)then
        if(xi(2)<3.5d0)then
            xi(2) = 3.5d0
        endif
        call TwoDPES(min(xi(2), 10.5d0), pii * xi(3) / 180d0, temp3)
        SO_corr = temp3
    endif
    !write(77,*) SO_corr

    if(xi(2)<5d0) then
        out_pot = out_pot + SO_corr
        if(out_pot>20990.11d0)then
            out_pot = 20990.11d0
        endif
        return
    endif

    call vpot(xi(1) / bohr, xi(2), xi(3), vx, vxso)

    if(SO_flag>0)then
        SO_corr2 = vxso - vx - 79.6d0
    endif
    vx = vx + 9354.59d0 + SO_corr2

    vx = vx - 2388.564169d0 * exp(-1d0 * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx + 18086.977116d0 * exp(-(1.307d0**1d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx - 71760.197585d0 * exp(-(1.307d0**2d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx + 154738.09175d0 * exp(-(1.307d0**3d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx - 215074.85646d0 * exp(-(1.307d0**4d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx + 214799.54567d0 * exp(-(1.307d0**5d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx - 148395.4285d0 * exp(-(1.307d0**6d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx + 73310.781453d0 * exp(-(1.307d0**7d0) * 0.785 * xi(1)**2) * 219.47463d0
    vx = vx + 42030.046d0

    SS = (1 - tanh(2.5d0 * (xi(2) - 8d0))) / 2d0
    !write(77,*) V,SO_corr
    out_pot = out_pot + SO_corr

    out_pot = SS * out_pot + (1d0 - SS) * vx

    if(out_pot>20990.11d0)then
        out_pot = 20990.11d0
        ! write(*,*) 'hit ceiling'
    endif

    contains
        ! we've moved the inner function inside, since this allows us to reuse parent-scope pieces

        function func_actual(xi)
            !    USE dynamic_parameters
            implicit none
            REAL(kind = 8), DIMENSION(:), INTENT(IN) :: xi
            REAL(kind = 8) :: func_actual
            integer :: i, j, k, ipp, jpp, ip, quitt, l1, l2, l3, l4, count2, l, jp, jj, kk, R, M
            integer :: count
            real(kind = 8) :: temp, weight, norm, somme, jac3(3), jac4(3), temp1, temp2, diff(3), pii, x, dx(3)
            real(kind = 8) :: ind7(count3)
            integer :: ind8(count3)

            pii = acos(-1d0)
            jac3 = xi

            call get_weights(jac3, count3, ind7, ind8)

!            ind7 = stored_weights
!            ind8 = stored_weights_ind

            quitt = 0
            norm = 0d0
            do ip = 1, count3
                quitt = quitt + 1
                norm = norm + ind7(ind8(count3 + 1 - ip))
                !   write(444,*) ind8(count3+1-ip),ind7(ind8(count3+1-ip))
                if(norm / ind7(ind8(count3 + 1 - ip))>1d8) goto 12
            enddo

            12 jac3 = xi
            Jac4 = jac3
            !norm=0d0
            temp = 0d0
            do i = 1, quitt
                jj = ind8(count3 + 1 - i)
                weight = ind7(ind8(count3 + 1 - i))
                !   norm=norm+weight
                count = 0
                do j = 0, order_1
                    do jpp = 0, order_1
                        do k = 0, order_1
                            if((j==0.and.jpp==0).or.(j==0.and.k==0).or.(jpp==0.and.k==0))then
                                count = count + 1
                                call basis_fnc(jac4, j, jpp, k, x, dx)
                                temp = temp + weight * b2(count, jj) * x
                            elseif(j==0.or.jpp==0.or.k==0)then
                                if(j + jpp<order_2 + 1.and.j + k<order_2 + 1.and.jpp + k<order_2 + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2(count, jj) * x
                                endif
                            else
                                if(j + jpp + k<order_3 + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2(count, jj) * x
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
            func_actual = temp / norm
!            deallocate(ind7, ind8)
            return
        end function func_actual

        function func_actual_min(xi)
            !    USE dynamic_parameters
            implicit none
            REAL(kind = 8), DIMENSION(:), INTENT(IN) :: xi
            REAL(kind = 8) :: func_actual_min
            integer :: i, j, k, ipp, jpp, ip, quitt, l1, l2, l3, l4, count2, l, jp, jj, kk, R, M
            integer :: count
            real(kind = 8) :: temp, weight, norm, somme, jac3(3), jac4(3), temp1, temp2, diff(3), pii, x, dx(3)
            real(kind = 8) :: ind7(count3)
            integer :: ind8(count3)
            pii = acos(-1d0)
            jac3 = xi

            call get_weights(jac3, count3, ind7, ind8)

            quitt = 0
            do ip = 1, count3
                if(ind7(ind8(count3)) / ind7(ind8(count3 + 1 - ip))>1d8) goto 12
                quitt = quitt + 1
            enddo
            !! outputs how many expansions are included in interpolation
            !write(701,*) quitt

            !!!!
            12 jac3 = xi
            Jac4 = jac3

            norm = 0d0
            temp = 0d0
            do i = 1, quitt
                jj = ind8(count3 + 1 - i)
                !   if(pot(jj)<E_limit)then
                weight = ind7(ind8(count3 + 1 - i))
                !   write(665,*) weight
                norm = norm + weight
                !   temp=temp+weight*b2_minimal(1,jj)
                count = 0
                do j = 0, order_1_min
                    do jpp = 0, order_1_min
                        do k = 0, order_1_min
                            if((j==0.and.jpp==0).or.(j==0.and.k==0).or.(jpp==0.and.k==0))then
                                count = count + 1
                                call basis_fnc(jac4, j, jpp, k, x, dx)
                                temp = temp + weight * b2_minimal(count, jj) * x
                            elseif(j==0.or.jpp==0.or.k==0)then
                                if(j + jpp<order_2_min + 1.and.j + k<order_2_min + 1.and.jpp + k<order_2_min + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2_minimal(count, jj) * x
                                endif
                            else
                                if(j + jpp + k<order_3_min + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2_minimal(count, jj) * x
                                endif
                            endif
                        enddo
                    enddo
                enddo

            enddo

            func_actual_min = temp / norm
!            deallocate(ind7, ind8)
            return
        end function func_actual_min

        function func_actual_seed(xi)
            !    USE dynamic_parameters
            implicit none
            REAL(kind = 8), DIMENSION(:), INTENT(IN) :: xi
            REAL(kind = 8) :: func_actual_seed
            integer :: i, j, k, ipp, jpp, ip, quitt, l1, l2, l3, l4, count2, l, jp, jj, kk, R, M
            integer :: count
            real(kind = 8) :: temp, weight, norm, somme, jac3(3), jac4(3), temp1, temp2, diff(3), pii, x, dx(3)
            real(kind = 8) :: ind7(count3)
            integer :: ind8(count3)
            pii = acos(-1d0)
            jac3 = xi
            !    allocate(ind7(count_seed),ind8(count_seed))
            count = 0
            do ip = 1, count_seed
                count = count + 1
                Jac4 = coords_seed(ip, :)
                call dist_metric(jac3, jac4, lab, somme)
                somme = somme**2
                ind7(count) = exp(-((somme) / d_seed(ip)**2)) / (((somme) / d_seed(ip)**2)**(zz_low / 2) + epss)
            enddo
            call indexxy(count_seed, ind7, ind8)

            quitt = 0
            norm = 0d0
            do ip = 1, count_seed
                quitt = quitt + 1
                norm = norm + ind7(ind8(count_seed + 1 - ip))
                if(norm / ind7(ind8(count_seed + 1 - ip))>1d8) goto 12
            enddo

            12 jac3 = xi
            Jac4 = jac3
            !norm=0d0
            temp = 0d0
            do i = 1, quitt
                jj = ind8(count_seed + 1 - i)
                weight = ind7(ind8(count_seed + 1 - i))
                !   norm=norm+weight
                count = 0
                do j = 0, order_1_min
                    do jpp = 0, order_1_min
                        do k = 0, order_1_min
                            if((j==0.and.jpp==0).or.(j==0.and.k==0).or.(jpp==0.and.k==0))then
                                count = count + 1
                                call basis_fnc(jac4, j, jpp, k, x, dx)
                                temp = temp + weight * b2_seed(count, jj) * x
                            elseif(j==0.or.jpp==0.or.k==0)then
                                if(j + jpp<order_2_min + 1.and.j + k<order_2_min + 1.and.jpp + k<order_2_min + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2_seed(count, jj) * x
                                endif
                            else
                                if(j + jpp + k<order_3_min + 1)then
                                    count = count + 1
                                    call basis_fnc(jac4, j, jpp, k, x, dx)
                                    temp = temp + weight * b2_seed(count, jj) * x
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo

            func_actual_seed = temp / norm
!            deallocate(ind7, ind8)
            return
        end function func_actual_seed

        subroutine get_weights(jac3, count3, stored_weights, stored_weights_ind)
            !  use dynamic_parameters
            implicit none

            integer :: i, ip, count
            real(kind = 8), intent(in) :: jac3(3)
            integer, intent(in) :: count3
            real(kind = 8), intent(inout) :: stored_weights(:)
            integer, intent(inout) :: stored_weights_ind(:)

            real(kind = 8) :: somme, jac4(3)
!            somme = 0d0
!            do i = 1, 3
!                somme = somme + (current_geom(i) - jac3(i))**2
!            enddo
!            somme = sqrt(somme)
!            if(somme<1d-11)then
!                return
!            else
!            current_geom = jac3
            count = 0
            do ip = 1, count3
                count = count + 1
                Jac4 = coords(ip, :)
                call dist_metric(jac3, jac4, lab, somme)
                somme = somme**2
                stored_weights(count) = exp(-((somme) / d(ip)**2)) / (((somme) / d(ip)**2)**(zz / 2) + epss)
            enddo
            call indexxy(count3, stored_weights, stored_weights_ind)
!            endif
            return
        end subroutine get_weights

        subroutine basis_fnc(jac4, l1, l2, l3, x, dx)
            !  use dynamic_parameters
            implicit none
            integer :: l1, l2, l3
            real(kind = 8) :: jac3(3), jac4(3), x, dx(3), PN(0:order_1), PD(0:order_1), temp1, temp2, a, b, p1, p2
            jac3 = jac4
            call LPN(order_1, jac3(3), PN, PD)
            x = (exp(alpha(1) * jac3(1))**l1) * (exp(alpha(2) * jac3(2))**l2) * PN(l3)
            dx(1) = dble(l1) * alpha(1) * (exp(alpha(1) * jac3(1))**l1) * (exp(alpha(2) * jac3(2))**l2) * PN(l3)
            dx(2) = dble(l2) * alpha(2) * (exp(alpha(1) * jac3(1))**l1) * (exp(alpha(2) * jac3(2))**l2) * PN(l3)
            dx(3) = (exp(alpha(1) * jac3(1))**l1) * (exp(alpha(2) * jac3(2))**l2) * PD(l3)
            return
        end subroutine basis_fnc

end subroutine IMLS_core

subroutine basis_size(D, order_1, order_2, order_3, lab, basis)
    implicit none

    integer, intent(in) :: D, order_1, order_2, order_3, lab
    integer, intent(out) :: basis
    integer :: count, ipp, jpp, i, j, kpp, k

    print *, order_1, order_2, order_3

    count = 0
    do i = 0, order_1
        do j = 0, order_1
            do k = 0, order_1
                if((i==0.and.j==0).or.(i==0.and.k==0).or.(j==0.and.k==0))then
                    count = count + 1
                elseif(i==0.or.j==0.or.k==0)then
                    if(i + j<order_2 + 1.and.i + k<order_2 + 1.and.j + k<order_2 + 1)then
                        count = count + 1
                    endif
                else
                    if(i + j + k<order_3 + 1)then
                        count = count + 1
                    endif
                endif
            enddo
        enddo
    enddo
    basis = count
    return
end subroutine basis_size

subroutine dist_metric(jac, jac2, lab, dist)
    integer :: i, lab
    real(kind = 8) :: jac(3), jac2(3), dist, temp(3), theta, theta2
    temp(1) = ((jac(1) - jac2(1)))**2
    temp(2) = ((jac(2) - jac2(2)))**2
    if(lab<3)then
        if(jac(3)>1d0)then
            jac(3) = 2d0 - jac(3)
        endif
        if(jac(3)<-1d0)then
            jac(3) = -2d0 - jac(3)
        endif
        if(jac2(3)>1d0)then
            jac2(3) = 2d0 - jac2(3)
        endif
        if(jac2(3)<-1d0)then
            jac2(3) = -2d0 - jac2(3)
        endif
        theta = acos(jac(3))
        theta2 = acos(jac2(3))
        !!$   temp(3)=((jac(3)-jac2(3)))**2*sqrt(jac(1)*jac2(1)*jac(2)*jac2(2))
        temp(3) = (theta - theta2)**2 * sqrt(jac(1) * jac2(1) * jac(2) * jac2(2))
    else
        temp(3) = ((jac(3) - jac2(3)))**2
    endif
    dist = 0d0
    do i = 1, 3
        dist = dist + temp(i)
    enddo
    dist = sqrt(dist)
    return
end subroutine dist_metric

subroutine INT_Cart(internal_temp, cart, mass, lab)
    implicit none
    integer :: lab
    real(kind = 8) :: internal(3), internal_temp(3), cart(9), mass(3), theta, cm, mu, M, d, hyper1, hyper2, hyper3
    !!! lab=1 (Jacobi, r,R,cos(theta)), lab=2 (valence, r1,r2,cos(theta)), lab=3
    !(bond dist, 12,13,23)
    !!now using theta instead of cos(theta)

    cart = 0d0
    internal = internal_temp

    if(lab==1)then
        if(internal(3)>1d0)then
            internal(3) = 2d0 - internal(3)
        endif
        if(internal(3)<-1d0)then
            internal(3) = -2d0 - internal(3)
        endif
        cm = mass(3) * internal(1) / (mass(2) + mass(3))
        theta = dacos(internal(3))
        cart(7) = internal(1)
        cart(1) = cm - internal(3) * internal(2)
        cart(2) = dsin(theta) * internal(2)
        !!$     cart(1)=cm-dcos(internal(3))*internal(2)
        !!$     cart(2)=dsin(internal(3))*internal(2)

    elseif(lab==2)then
        if(internal(3)>1d0)then
            internal(3) = 2d0 - internal(3)
        endif
        if(internal(3)<-1d0)then
            internal(3) = -2d0 - internal(3)
        endif
        cart(1) = internal(1)
        theta = dacos(internal(3))
        cart(7) = internal(3) * internal(2)
        cart(8) = dsin(theta) * internal(2)
        !!$     cart(7)=dcos(internal(3))*internal(2)
        !!$     cart(8)=dsin(internal(3))*internal(2)

    elseif(lab==3)then
        if(internal(1)>internal(2) + internal(3).or.internal(2)>internal(3) + internal(1).or.internal(3)>internal(1) + internal(2))then
            write(*, *) internal(:), 'triangle condition not met'
        endif
        cart(4) = internal(1)
        cart(7) = (internal(1)**2 + internal(2)**2 - internal(3)**2) / (2d0 * internal(1))
        cart(8) = sqrt(internal(2)**2 - cart(7)**2)
        if(abs(cart(7))>internal(2))then
            cart(8) = 0d0
        endif
    else
        M = mass(1) + mass(2) + mass(3)
        mu = sqrt(mass(1) * mass(2) * mass(3) / M)

        d = sqrt((mass(1) / mu) * (1d0 - mass(1) / M))
        hyper1 = (internal(1) / sqrt(2d0)) * sqrt(1d0 + dsin(internal(2)) * dcos(internal(3)))
        hyper2 = (internal(1) / sqrt(2d0)) * sqrt(1d0 - dsin(internal(2)) * dcos(internal(3)))
        hyper3 = dsin(internal(2)) * dsin(internal(3)) / sqrt(1d0 - dsin(internal(2))**2 * dcos(internal(3))**2)

        internal(1) = d * hyper2
        internal(2) = (1d0 / d) * hyper1
        internal(3) = dacos(hyper3)
        cm = mass(3) * internal(1) / (mass(2) + mass(3))
        !     theta=dacos(internal(3))

        cart(7) = internal(1)
        cart(1) = cm - dcos(internal(3)) * internal(2)!cm-internal(3)*internal(2)
        cart(2) = dsin(internal(3)) * internal(2)
    endif
    return
end subroutine INT_Cart

subroutine Cart_INT(internal, cart, mass, lab)

    implicit none
    integer :: i, lab
    real(kind = 8) :: internal(3), cart(9), mass(3), theta, cm(3), dist, vec(3), vec2(3), mu, M, d, hyper1, hyper2, hyper3, x, y
    !!!using theta, not cos(theta)
    internal = 0d0
    if(lab==1)then
        internal(1) = sqrt((cart(4) - cart(7))**2 + (cart(5) - cart(8))**2 + (cart(6) - cart(9))**2)
        dist = 0d0
        do i = 1, 3
            vec(i) = cart(3 + i) - cart(6 + i)
            cm(i) = (cart(3 + i) * mass(2) + cart(6 + i) * mass(3)) / (mass(2) + mass(3))
            vec2(i) = cart(i) - cm(i)
            dist = dist + vec(i) * vec2(i)
        enddo
        internal(2) = sqrt((cart(1) - cm(1))**2 + (cart(2) - cm(2))**2 + (cart(3) - cm(3))**2)
        internal(3) = dist / (internal(1) * internal(2))
        !!$     internal(3)=dacos(dist/(internal(1)*internal(2)))

    elseif(lab==2)then
        internal(1) = sqrt((cart(1) - cart(4))**2 + (cart(2) - cart(5))**2 + (cart(3) - cart(6))**2)
        internal(2) = sqrt((cart(7) - cart(4))**2 + (cart(8) - cart(5))**2 + (cart(9) - cart(6))**2)
        dist = 0d0
        do i = 1, 3
            vec(i) = cart(i) - cart(3 + i)
            vec2(i) = cart(6 + i) - cart(3 + i)
            dist = dist + vec(i) * vec2(i)
        enddo
        internal(3) = dist / (internal(1) * internal(2))
        !!$     internal(3)=dacos(dist/(internal(1)*internal(2)))
    elseif(lab==3)then
        internal(1) = sqrt((cart(1) - cart(4))**2 + (cart(2) - cart(5))**2 + (cart(3) - cart(6))**2)
        internal(2) = sqrt((cart(1) - cart(7))**2 + (cart(2) - cart(8))**2 + (cart(3) - cart(9))**2)
        internal(3) = sqrt((cart(7) - cart(4))**2 + (cart(8) - cart(5))**2 + (cart(9) - cart(6))**2)
    else
        !!!!hyperspherical
        !!! start with scaled Jacobi
        internal(1) = sqrt((cart(4) - cart(7))**2 + (cart(5) - cart(8))**2 + (cart(6) - cart(9))**2)
        dist = 0d0
        do i = 1, 3
            vec(i) = cart(3 + i) - cart(6 + i)
            cm(i) = (cart(3 + i) * mass(2) + cart(6 + i) * mass(3)) / (mass(2) + mass(3))
            vec2(i) = cart(i) - cm(i)
            dist = dist + vec(i) * vec2(i)
        enddo
        internal(2) = sqrt((cart(1) - cm(1))**2 + (cart(2) - cm(2))**2 + (cart(3) - cm(3))**2)

        M = mass(1) + mass(2) + mass(3)
        mu = sqrt(mass(1) * mass(2) * mass(3) / M)
        d = sqrt((mass(1) / mu) * (1d0 - mass(1) / M))
        internal(1) = (1d0 / d) * internal(1)
        internal(2) = d * internal(2)
        internal(3) = dacos(dist / (internal(1) * internal(2)))
        !! start hyper
        hyper1 = sqrt(internal(1)**2 + internal(2)**2)
        hyper2 = datan(sqrt((internal(2)**2 - internal(1)**2)**2 + (2d0 * dist)**2) / (2d0 * internal(1) * internal(2) * dsin(internal(3))))
        y = 2d0 * dist / sqrt((internal(2)**2 - internal(1)**2)**2 + (2d0 * dist)**2)
        x = (internal(2)**2 - internal(1)**2) / sqrt((internal(2)**2 - internal(1)**2)**2 + (2d0 * dist)**2)
        hyper3 = datan2(y, x)
        internal(1) = hyper1
        internal(2) = hyper2
        internal(3) = hyper3
    endif
    return
end subroutine Cart_INT

subroutine LPN(N, X, PN, PD)
    implicit none
    integer :: N, K
    real(kind = 8) :: X, PN(0:N), PD(0:N), P0, P1, PF
    PN(0) = 1.0D0
    PN(1) = X
    PD(0) = 0.0D0
    PD(1) = 1.0D0
    P0 = 1.0D0
    P1 = X
    DO K = 2, N
        PF = (2.0D0 * dble(K) - 1.0D0) / dble(K) * X * P1 - (dble(K) - 1.0D0) / dble(K) * P0
        PN(K) = PF
        IF (DABS(X)==1.0D0) THEN
            PD(K) = 0.5D0 * X**dble(K + 1) * dble(K) * (dble(K) + 1.0D0)
        ELSE
            PD(K) = dble(K) * (P1 - X * PF) / (1.0D0 - X * X)
        ENDIF
        P0 = P1
        P1 = PF
    ENDDO
end SUBROUTINE LPN
