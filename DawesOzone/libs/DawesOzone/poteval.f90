module poteval
    implicit none

    integer :: Nd, Ld, nfj, lmax, Nip
    parameter(Nd = 9, Ld = 9, nfj = 9, lmax = 2, Nip = 1)

contains

    subroutine cnjel(roo, ic, ipp, cnjj)
        use splin

        implicit none

        !        parameter (lmax = 2, nfj = 9)
        integer :: ipp, ni, nj, nt, nf, nft, ms, ml, mj, mps, mpla, mpj, mma, &
                mmb, mm, m2, m1, lb, l, jt, jjp, jj, j, it, imm, is, i, il, ima, mla, &
                ic

        real(kind = 8) :: xs, xms, xml, xmj, xl, xj, xcn, xcleb, xcj, xci, roo
        real(kind = 8) :: clb(0:lmax), cn(0:lmax, 0:2 * lmax + 1), Cnmmp(4, 0:3, 0:3, 0:2), ci(nfj, 3)
        real(kind = 8) :: cnjj(nfj, nfj, 0:lmax, 0:2 * lmax + 1)

        integer :: ims(nfj, 3), iml(nfj, 3), jnt(nfj), ijj(nfj), imj(nfj)

        double precision :: coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
        integer :: ll(splin_nlt, 6)

        ! element de matrice de O2
        !      write(9,*)'ipp=',ipp
        if(ipp==1)then
            mmb = 0
        else
            print *, 'erreur de ipp', ipp
            stop
        endif
        ! Decomposition de la fonction d'onde |JMj> de l'atome en |LSMlMs>
        l = 1
        is = 1
        xl = dfloat(l)
        xs = dfloat(is)
        !     print *,'J Mj Ml Ms Clebsh'
        nf = 0
        do jj = iabs(l - is), l + is
            xj = dfloat(jj)
            do mj = -jj, jj
                nf = nf + 1
                ijj(nf) = jj
                imj(nf) = mj
                xmj = dfloat(mj)
                nt = 0
                do ml = -l, l
                    nt = nt + 1
                    xml = dfloat(ml)
                    ms = mj - ml
                    xms = dfloat(ms)
                    call cleb(xl, xs, xj, xml, xms, xmj, 1, xcleb)
                    !     print '(4f6.2,f16.8)',xj,xmj,xml,xms,xcleb
                    ci(nf, nt) = xcleb
                    ims(nf, nt) = ms
                    iml(nf, nt) = ml
                enddo
                jnt(nf) = nt
                !     print *
            enddo
        enddo
        nft = nf
        !     print *,'nft=',nft
        if(nft/=nfj)then
            print *, 'erreur de dimension de ci,ims et iml', nft
            stop
        endif
        ! Calcul des Cnelec(mla,mpla,ipp)
        do mla = -l, l
            do mpla = -l, l
                mma = mla - mpla
                call cnelr(coef, ll, roo, ic, mla, mpla, ipp, mma, mmb, clb)
                do lb = 0, lmax
                    m1 = mla + l
                    m2 = mpla + l
                    cnmmp(ipp, m1, m2, lb) = clb(lb)
                enddo
            enddo
        enddo
        ! Calcul des Cn pour les elements de matrice <JMj| |>< | |J'Mj'>
        do i = 1, nft
            jj = ijj(i)
            mj = imj(i)
            do j = 1, nft
                jjp = ijj(j)
                mpj = imj(j)
                !      write(9,100)jj,mj,jjp,mpj
                ni = nft - i + 1
                nj = nft - j + 1
                do il = 0, lmax
                    do mm = -il, il
                        imm = mm + il
                        if(imm>5)print *, 'erreur de dim de C6', imm
                        cn(il, imm) = 0.d0
                    enddo
                enddo
                do it = 1, jnt(i)
                    ms = ims(i, it)
                    mla = iml(i, it)
                    m1 = mla + l
                    xci = ci(i, it)
                    if(dabs(xci)>1.e-5)then
                        do jt = 1, jnt(j)
                            mps = ims(j, jt)
                            mpla = iml(j, jt)
                            m2 = mpla + l
                            xcj = ci(j, jt)
                            if(dabs(xcj)>1.e-5)then
                                if(iabs(ms - mps)<1.e-4)then
                                    !      write(9,99),ms,mla,mps,mpla
                                    do lb = 0, lmax
                                        ima = mla - mpla + lb
                                        xcn = cnmmp(ipp, m1, m2, lb)
                                        if(ima<0.and.dabs(xcn)>1.e-5)then
                                            print *, '|Ma| >  Lb : STOP', mla, mpla, lb, xcn
                                            stop
                                        endif
                                        !     print *,ipp,m1,m2,lb,xcn
                                        !     xgw,  fix cn(:,-1) bug
                                        if (ima >= 0) then
                                            cn(lb, ima) = cn(lb, ima) + xcn * xci * xcj
                                        end if
                                    enddo
                                endif
                            endif
                        enddo
                    endif
                enddo
                do il = 0, lmax
                    do mm = -il, il
                        imm = mm + il
                        if(dabs(cn(il, imm))>1.e-5)then
                            !      	write(9,101)il,mm,mmb,cn(il,imm)
                        endif
                        cnjj(ni, nj, il, imm) = cn(il, imm)
                    enddo
                enddo
            enddo
        enddo
        !    99  format('    Ms=', i3, ' Ml=', i3, ' Mps=', i3, ' Mpl=', i3)
        !    100  format('<J=', i3, ' Mj=', i3, ' Jp=', i3, ' Mjp=', i3, '>')
        !    101  format('Cnel(Lb=', i3, 'Ma=', i3, 'Mb=', i3, ')=', f16.8)
        return
    end subroutine cnjel

    subroutine c6r(roo, ic, mla, mpla, Ma, ipp, Mb, c6lb)
        use splin
        implicit none

        ! a partir du fichier tabCn.res
        !        implicit real(kind = 8)(a-h, o-z)
        !        parameter (Nip = 1, lmax = 2)
        real(kind = 8) :: roo, C6lb(0:lmax), xi(splin_nr), yi(splin_nr), b(splin_nr), c(splin_nr), d(splin_nr)
        real(kind = 8) :: z0, un, xl1, xl1p, xl2, xl2p, c6t

        integer :: ipp, l1, l1p, l2, l2p, n, lb, Ma, Mb, mla, ml, mpla, &
                ic, iMa, iMb, iLb, imla, impla, jpp, la, lambda

        double precision :: coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
        integer :: ll(splin_nlt, 6)

        l1 = 1
        l1p = 1
        l2 = 1
        l2p = 1

        z0 = 0.d0
        un = 1.d0

        !        common/droo/roo
        !        common/dcc/ic
        !        data l1, l1p/1, 1/
        !        data l2, l2p/1, 1/
        !        data z0, un/0.d0, 1.d0/
        n = l1 + l1p + l2 + l2p + 2
        xl1 = dfloat(l1)
        xl1p = dfloat(l1p)
        xl2 = dfloat(l2)
        xl2p = dfloat(l2p)
        !      write(6,101)l1,l1p,l2,l2p
        !        101   format('# l1 l1p l2 l2p =', 4i3)
        do lb = 0, lmax
            c6lb(lb) = 0.d0
        enddo
        do lb = iabs(Mb), lmax
            c6t = 0.d0
            do la = iabs(Ma), lmax
                ml = min0(la, lb)
                do lambda = iabs(la - lb), la + lb
                enddo
            enddo
            ic = ic + 1
            call initspline(ic, xi, yi, b, c, d, coef)
            c6t = ispline(roo, xi, yi, b, c, d, splin_nr)
            !      read(10,*)imla,impla,jpp,iMa,iMb,iLb,c6t
            imla = ll(ic, 1)
            impla = ll(ic, 2)
            jpp = ll(ic, 3)
            iMa = ll(ic, 4)
            iMb = ll(ic, 5)
            iLb = ll(ic, 6)
            !       write(11,106),imla,impla,jpp,iMa,iMb,iLb,c6t
            if(iabs(imla - mla)/=0)then
                print *, '# erreur de lecture de Mla', imla, mla
                stop
            endif
            if(iabs(impla - mpla)/=0)then
                print *, '# erreur de lecture de Mpla', impla, mpla
                stop
            endif
            if(iabs(ima - ma)/=0)then
                print *, '# erreur de lecture de Ma', ima, ma
                stop
            endif
            if(iabs(imb - mb)/=0)then
                print *, '# erreur de lecture de Mb', imb, mb
                stop
            endif
            if(iabs(ilb - lb)/=0)then
                print *, '# erreur de lecture de Lb', ilb, lb
                stop
            endif
            if(iabs(ipp - jpp)/=0)then
                print *, '# erreur de lecture de ipp', jpp, ipp
                stop
            endif
            c6lb(lb) = c6t
        enddo
        !        106  format(i3, '&',i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',g16.8)
        return
    end subroutine c6r

    subroutine c6jmat(roo, ic, ipp, c6jj)
        use splin

        !        implicit real(kind = 8)(a-h, o-z)
        !        parameter (lmax = 2, nfj = 9)
        implicit none

        integer :: ic, ipp, i, nft, il, ima, imm, is, it, j, jj, jjp, jt, l, lb, m1, m2, &
                mj, ml, mla, mm, mma, mmb, mpj, mpla, mps, ms, nf, ni, nj, nt
        real(kind = 8) :: xc6, xci, xcj, xcleb, xj, xl, xmj, xml, xms, xs, roo

        real(kind = 8) :: c6lb(0:lmax), c6(0:lmax, 0:2 * lmax + 1), ci(nfj, 3), C6mmp(4, 0:3, 0:3, 0:2)
        integer :: jnt(nfj), ijj(nfj), imj(nfj), ims(nfj, 3), iml(nfj, 3)
        real(kind = 8) :: c6jj(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        ! element de matrice de O2
        !      write(9,*)'ipp=',ipp
        if(ipp==1)then
            mmb = 0
        else
            print *, 'erreur de ipp', ipp
            stop
        endif
        ! Decomposition de la fonction d'onde |JMj> de l'atome  en |LSMlMs>
        l = 1
        is = 1
        xl = dfloat(l)
        xs = dfloat(is)
        !     print *,'J Mj Ml Ms Clebsh'
        nf = 0
        do jj = iabs(l - is), l + is
            xj = dfloat(jj)
            do mj = -jj, jj
                nf = nf + 1
                ijj(nf) = jj
                imj(nf) = mj
                xmj = dfloat(mj)
                nt = 0
                do ml = -l, l
                    nt = nt + 1
                    xml = dfloat(ml)
                    ms = mj - ml
                    xms = dfloat(ms)
                    call cleb(xl, xs, xj, xml, xms, xmj, 1, xcleb)
                    !     print '(4f6.2,f16.8)',xj,xmj,xml,xms,xcleb
                    ci(nf, nt) = xcleb
                    ims(nf, nt) = ms
                    iml(nf, nt) = ml
                enddo
                jnt(nf) = nt
                !     print *
            enddo
        enddo
        nft = nf
        !     print *,'nft=',nft
        if(nft/=nfj)then
            print *, 'erreur de dimension de ci,ims et iml', nft
            stop
        endif
        ! Calcul des C6(mla,mpla,ipp)
        do mla = -l, l
            do mpla = -l, l
                mma = mla - mpla
                call C6r(roo, ic, mla, mpla, mma, ipp, mmb, c6lb)
                do lb = 0, lmax
                    m1 = mla + l
                    m2 = mpla + l
                    c6mmp(ipp, m1, m2, lb) = c6lb(lb)
                enddo
            enddo
        enddo
        ! Calcul des C6 pour les elements de matrice <JMj| |>< | |J'Mj'>
        do i = 1, nft
            jj = ijj(i)
            mj = imj(i)
            do j = 1, nft
                jjp = ijj(j)
                mpj = imj(j)
                !      write(9,100)jj,mj,jjp,mpj
                ni = nft - i + 1
                nj = nft - j + 1
                do il = 0, lmax
                    do mm = -il, il
                        imm = mm + il
                        if(imm>5)print *, 'erreur de dim de C6', imm
                        c6(il, imm) = 0.d0
                    enddo
                enddo
                do it = 1, jnt(i)
                    ms = ims(i, it)
                    mla = iml(i, it)
                    m1 = mla + l
                    xci = ci(i, it)
                    if(dabs(xci)>1.e-5)then
                        do jt = 1, jnt(j)
                            mps = ims(j, jt)
                            mpla = iml(j, jt)
                            m2 = mpla + l
                            xcj = ci(j, jt)
                            if(dabs(xcj)>1.e-5)then
                                if(iabs(ms - mps)<1.e-4)then
                                    !      write(9,99),ms,mla,mps,mpla
                                    do lb = 0, lmax
                                        ima = mla - mpla + lb
                                        xc6 = c6mmp(ipp, m1, m2, lb)
                                        if(ima<0.and.dabs(xc6)>1.e-5)then
                                            print *, '|Ma| >  Lb : STOP', mla, mpla, lb, xc6
                                            stop
                                        endif
                                        !     print *,ipp,m1,m2,lb,xc6
                                        !     xgw,  fix c6(:,-1) bug
                                        if (ima >= 0) then
                                            c6(lb, ima) = c6(lb, ima) + xc6 * xci * xcj
                                        end if
                                    enddo
                                endif
                            endif
                        enddo
                    endif
                enddo
                do il = 0, lmax
                    do mm = -il, il
                        imm = mm + il
                        if(dabs(c6(il, imm))>1.e-5)then
                            !      	write(9,101)il,mm,mmb,c6(il,imm)
                        endif
                        c6jj(ni, nj, il, imm) = c6(il, imm)
                    enddo
                enddo
            enddo
        enddo
        !        99  format('Ms=', i3, ' Ml=', i3, ' Mps=', i3, ' Mpl=', i3)
        !        100  format('<J=', i3, ' Mj=', i3, ' Jp=', i3, ' Mjp=', i3, '>')
        !        101  format('C6(Lb=', i3, 'Ma=', i3, 'Mb=', i3, ')=', f16.8)
        return
    end subroutine c6jmat

    subroutine cleb(a, b, c, d, e, f, kselc, rac)
        !#######################################################################
        !     taro tamura, c.p.c. 1(1970)337-342.  (correcte4 version 19/1171@ #
        !     calculates clebsch-gordon coefficients.                          #
        !     see g.racah, phys.rev. 62(1942)438.                              #
        ! --- this subroutine has been modified according to the comments made #
        !     by j. g. wills in c.p.c. 2(1971)381-382.                         #
        !----------------------------------------------------------------------#
        !     rac = clebsh (j1=a, j2=b, j3=c; m1=d, m2=e, m3=f)                #
        !     kselc # 0  : test for selection rule before calculation          #
        !#######################################################################
        use splin

        implicit none
        !        double precision (a-h, o-z)

        !        parameter (nfctmx = 250)
        !        common/logfac/facl(nfctmx)
        !        data ntimes/0/

        integer :: nfctmx, ntimes, i, ia, ib, ic, id, ie, if, ig4, kselc, &
                i1, i2, i3, i4, i5, i6, i7, i8, iabc, iabcp, iamd, iapd, ibca, &
                ibme, ibpe, icab, icmf, icpf, ig, ig2, j, nterm, nz, nzm1, nzmi, &
                nzmic2, nzmic3, nzmx, nzt1, nzt2, nzt3, nzt4, nzt5
        double precision :: facl(250), rac, f1, a, b, c, d, e, f, ai, a1, &
                f2, fb, fc2, s1, sqfclg, termlg, x, y

        nfctmx = 250
        ntimes = 0

        if (ntimes ==0) call faclog(facl, nfctmx)
        ntimes = ntimes + 1
        !
        ia = ifix(sngl(2. * a))
        ib = ifix(sngl(2. * b))
        ic = ifix(sngl(2. * c))
        id = ifix(sngl(2. * d))
        ie = ifix(sngl(2. * e))
        if = ifix(sngl(2. * f))
        !
        rac = 0.0
        ig4 = ia + ib + ic
        if(kselc==0) go to 170
        if(id + ie/=if) go to 1000
        if(mod(ig4, 2)/=0) go to 1000
        if(ia + ib - ic<0.or.ic - iabs(ia - ib)<0) go to 1000
        if(min0(ia - iabs(id), ib - iabs(ie), ic - iabs(if))<0) go to 1000
        if(mod(ib + ie, 2)/=0.or.mod(ic + if, 2)/=0) go to 1000
        170 if(ia==0.or.ib==0) go to 175
        if(ic < 0) then
            go to 1000
        elseif (ic == 0) then
            go to 180
        else
            go to 185
        end if
        175 rac = 1.0
        go to 1000
        180 fb = ib + 1
        rac = ((-1.0)**((ia - id) / 2)) / dsqrt(fb)
        go to 1000
        185 if(id/=0.or.ie/=0) go to 250
        ig2 = ig4 / 2
        if(mod(ig2, 2)/=0) go to 1000
        ig = ig2 / 2
        i1 = (ia + ib - ic) / 2 + 1
        i2 = (ic + ia - ib) / 2 + 1
        i3 = (ic + ib - ia) / 2 + 1
        i4 = ig2 + 2
        i5 = ig + 1
        i6 = (ig2 - ia) / 2 + 1
        i7 = (ig2 - ib) / 2 + 1
        i8 = (ig2 - ic) / 2 + 1
        f1 = dexp(.5 * (facl(i1) + facl(i2) + facl(i3) - facl(i4))&
                + (facl(i5) - facl(i6) - facl(i7) - facl(i8)))
        f2 = ic + 1
        f2 = dsqrt(f2)
        s1 = 1 - 2 * mod((ig2 + ic) / 2, 2)
        rac = s1 * f1 * f2
        go to 1000
        250 fc2 = ic + 1
        iabcp = ig4 / 2 + 1
        iabc = iabcp - ic
        icab = iabcp - ib
        ibca = iabcp - ia
        iapd = (ia + id) / 2 + 1
        iamd = iapd - id
        ibpe = (ib + ie) / 2 + 1
        ibme = ibpe - ie
        icpf = (ic + if) / 2 + 1
        icmf = icpf - if
        sqfclg = 0.5 * (dlog(fc2) - facl(iabcp + 1)&
                + facl(iabc) + facl(icab) + facl(ibca)&
                + facl(iapd) + facl(iamd) + facl(ibpe)&
                + facl(ibme) + facl(icpf) + facl(icmf))
        nzmic2 = (ib - ic - id) / 2
        nzmic3 = (ia - ic + ie) / 2
        nzmi = max0(0, nzmic2, nzmic3) + 1
        nzmx = min0(iabc, iamd, ibpe)
        s1 = (-1.0)**(nzmi - 1)
        nz = nzmi
        nzm1 = nzmi - 1
        nzt1 = iabc - nzm1
        nzt2 = iamd - nzm1
        nzt3 = ibpe - nzm1
        nzt4 = nz - nzmic2
        nzt5 = nz - nzmic3
        termlg = sqfclg - facl(nz) - facl(nzt1) - facl(nzt2)&
                - facl(nzt3) - facl(nzt4) - facl(nzt5)
        rac = s1 * dexp(termlg)
        !
        nterm = nzmx - nzmi + 1
        a1 = 1.
        i = nterm
        do j = 1, nterm
            x = float((iabc - nzm1 - i) * (iamd - nzm1 - i) * (ibpe - nzm1 - i))
            y = float((nzm1 + i) * (-nzmic2 + nzm1 + i) * (-nzmic3 + nzm1 + i))
            a1 = 1. - a1 * x / y
            i = i - 1
        end do
        rac = rac * a1
        if(dabs(rac)<0.000001) rac = 0.0
        !
        1000 return
    end

    subroutine faclog(fct, nfctmx)
        !#######################################################################
        !#    initialisation of logarithms of factorials array                 #
        !#######################################################################
        use splin

        implicit none

        integer :: nfctmx, ntimes, i
        double precision :: fct(nfctmx), ai
        !        parameter (nfctmx = 1001)
        !        common /logfac/ fct(nfctmx)
        !        data ntimes /0/
        !
        ntimes = ntimes + 1
        if (ntimes > 1) return
        fct(1) = 0.d0
        do i = 1, nfctmx - 1
            ai = i
            fct(i + 1) = fct(i) + dlog(ai)
        end do
        !
        return
    end

    subroutine cnelr(coef, ll, roo, ic, mla, mpla, ipp, Ma, Mb, clb)
        use splin
        ! energie electrostatique de X(3P)+OH(X2Pi)
        implicit real(kind=8)(a-h, o-z)

        integer :: ic, Ma, Mb, mla, mpla, lb, l1, l2, ipp, jpp, &
                iLb, iMa, iMb, imla, impla

!        parameter (Nip = 1, lmax = 2)
        dimension clb(0:lmax)
        dimension xi(splin_nr), yi(splin_nr), b(splin_nr), c(splin_nr), d(splin_nr)

        double precision :: coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
        integer :: ll(splin_nlt, 6)

!        common/droo/roo
!        common/dcc/ic
        if(Ma/=mla - mpla)then
            write(6, *)'Erreur de Ma', Ma
            stop
        endif
        do lb = 1, lmax
            clb(lb) = 0.d0
        enddo
        l1 = 2
        do l2 = 1, lmax
            Lb = l2
            if(iabs(Mb)<=Lb)then
                ic = ic + 1
                call initspline(ic, xi, yi, b, c, d, coef)
                cc = ispline(roo, xi, yi, b, c, d, splin_nr)
                !      read(10,108)imla,impla,jpp,iMa,iMb,iLb,cc
                imla = ll(ic, 1)
                impla = ll(ic, 2)
                jpp = ll(ic, 3)
                iMa = ll(ic, 4)
                iMb = ll(ic, 5)
                iLb = ll(ic, 6)
                !       write(11,108) imla,impla,jpp,iMa,iMb,iLb,cc
                if(iabs(imla - mla)/=0)then
                    print *, '# erreur de lecture de Mla', imla, mla
                    stop
                endif
                if(iabs(impla - mpla)/=0)then
                    print *, '# erreur de lecture de Mpla', impla, mpla
                    stop
                endif
                if(iabs(ima - ma)/=0)then
                    print *, '# erreur de lecture de Ma', ima, ma
                    stop
                endif
                if(iabs(imb - mb)/=0)then
                    print *, '# erreur de lecture de Mb', imb, mb
                    stop
                endif
                if(iabs(ilb - lb)/=0)then
                    print *, '# erreur de lecture de Lb', ilb, lb
                    stop
                endif
                if(iabs(ipp - jpp)/=0)then
                    print *, '# erreur de lecture de ipp', jpp, ipp
                    stop
                endif
!                108  format(i3, '&',i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',g16.8)
                !      print *
                clb(l2) = cc
            else
                clb(l2) = 0.d0
            endif
        enddo
        return
    end

    !       program VLGRoo
    !       use splin
    !       implicit real*8(a-h,o-z)
    ! calcul de l'energie  Electrostatique +  Dispersion pour l'etat fondamental
    ! gam en degree, r et R en bohr, energies en cm-1
    ! pour 1.8 < r < 3.3 bohr
    !       data tocm/219474.63d0/
    !       common/drob/q20O2,rooi
    !       pi=dacos(-1.d0)
    !       torad=pi/180.d0
    !       rooi=0.d0
    ! Lecture des data
    !       call inpdat
    !
    !       write(6,*)'# r=dist(O-O)   R=dist(O-O2), gam=angle(OG) '
    !       write(6,*)'# roo(a0)  R(a0)  Gam(deg)  En(X)   En(X+SO) (cm-1)'
    !       do roo = 1.8,3.3,0.2
    !       do gam=0.d0,90.d0,30.d0
    !       do rr=8.d0,20.d0,0.2d0
    !       call vpot(roo,rr,gam,vx,vxso)
    !      write(7,*)
    !       write(7,*)'Apres Diagonalisation'
    !       write(6,'(5g16.8)')roo,rr,gam,vx,vxso
    !       enddo
    !       write(6,*)
    !       enddo
    !       enddo
    !       end
    subroutine buildmatrix(&
            roo, rr, gam, vmf, vmfso, &
            cnee1, c6jj1)
        use splin
        implicit none
        ! 7 constantes: q_1^0(OH), q_2^0(OH),q_2^2(OH),Q_0(O),DelE(3P_0-^3P_1)
        !  Del1O=DE(^3P_1-^3P_0), Del2O=DE(^3P_2-^3P_0),DelOH=DE(^2Pi_1/2-^2Pi_3/2)

        real(kind = 8) :: roo, rr, gam

        real(kind = 8) :: VM(Nd, Nd), VMF(Nd, Nd), VMFSO(Nd, Nd)
        real(kind = 8) :: E(9, 9), F(9, 9), G(9, 9)
        real(kind = 8) :: Q20(9, 9), Q21(9, 9), Q22(9, 9), Q2m1(9, 9), Q2m2(9, 9)
        real(kind = 8) :: A(9, 9), B(9, 9), C(9, 9)
        real(kind = 8) :: EFS(Nd, Nd), E1(Nd, Nd), E2(Nd, Nd)
        real(kind = 8) :: xi(splin_nr), yi(splin_nr), bq(splin_nr), cq(splin_nr), dq(splin_nr)
        real(kind = 8) :: cnee1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        real(kind = 8) :: c6jj1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)

        real(kind = 8) :: a1, a2, autocm, autoev, c2, cc, dd, del1O, del2O, gm, &
                pi, Q0O, q10O2, q20O2, rooi, s, s2, torad, xc6, xcn, xe, xe2
        integer :: i, il, imb, nfj, imm, imm2, ipp, j, mm, mmb

        double precision :: coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
        integer :: ll(splin_nlt, 6)


        !        common/cne/cnee1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        !        common/c6j/c6jj1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        !        common/drob/q20O2, rooi
        ! Moments en a.u. ; q20_OH=Q20(OH)
        Q0O = -1.90d0
        !      data q10O2,q20O2/0.0d0,-0.4546d0/
        q10O2 = 0.0d0
        ! Energies en cm-1
        del1O = 158.265d0
        del2O = 226.977d0
        autocm = 219474.625d0

        !      ! data autocm/2.194746313708d+5/
        ! alpha en a0^3, energie d'ionisation in eV
        autoev = 27.2113957d0
        !      ! autoev = 27.21138505d0
        pi = dacos(-1.d0)
        torad = pi / 180.d0
        gm = gam * torad
        ! Initialisation des Cn et q20 en fct de roo
        !        if(dabs(roo - rooi)>1.e-06) then
        call init(roo, cnee1, c6jj1, coef, ll)
        ! evaluation de q20O2  pour roo
        call initspline2(xi, yi, bq, cq, dq, rro)
        q20O2 = ispline(roo, xi, yi, bq, cq, dq, splin_nr)
        !       write(11,*)roo,q10O2,q20O2
        !       write(11,*)
        !        endif
        rooi = roo
        ! construction de la matrice
        ! element <+1| |+1>
        ipp = 1
        mmb = 0
        !      print *,'ipp=',ipp
        do i = 1, nfj
            do j = 1, nfj
                xe = 0.d0
                xe2 = 0.d0
                do il = iabs(mmb), lmax
                    do mm = -il, il
                        imm = mm + il
                        imm2 = -mm + il
                        imb = mmb + il
                        call matd(il, imm2, imb, gm, dd)
                        xcn = cnee1(i, j, il, imm)
                        xe = xe + xcn * dd / rr**(il + 3)
                        xc6 = c6jj1(i, j, il, imm)
                        xe2 = xe2 + xc6 * dd / rr**6
                        !      if(dabs(xcn).gt.1.e-5)then
                        !      print '(4i3)',i,j,il,mm
                        !      print '(4f16.8)',gm,xcn,dd,xe
                        !      endif
                    enddo
                enddo
                E1(i, j) = xe
                E2(i, j) = xe2
            enddo
        enddo
        ! Construction des matrices A,B,C,E,F,G,EFS
        s = dsin(gm)
        cc = dcos(gm)
        s2 = dsin(gm / 2.d0)
        c2 = dcos(gm / 2.d0)
        do i = 1, 9
            do j = 1, 9
                E(i, j) = 0.d0
                F(i, j) = 0.d0
                G(i, j) = 0.d0
                EFS(i, j) = 0.d0
            enddo
        enddo
        do i = 1, Nd
            do j = 1, Nd
                vm(i, j) = 0.d0
            enddo
        enddo
        ! Matrix EFS (Structure Fine)
        EFS(1, 1) = 0.d0
        EFS(2, 2) = 0.d0
        EFS(3, 3) = 0.d0
        EFS(4, 4) = 0.d0
        EFS(5, 5) = 0.d0
        EFS(6, 6) = del1O
        EFS(7, 7) = del1O
        EFS(8, 8) = del1O
        EFS(9, 9) = del2O
        ! Matrix E
        E(1, 1) = -1.d0
        E(2, 2) = 0.5d0
        E(3, 3) = 1.d0
        E(4, 4) = 0.5d0
        E(5, 5) = -1.d0
        E(6, 6) = 0.5d0
        E(7, 7) = -1.d0
        E(8, 8) = 0.5d0
        E(2, 6) = -1.5d0
        E(3, 9) = -dsqrt(2.d0)
        E(4, 8) = 1.5d0
        E(6, 2) = -1.5d0
        E(8, 4) = 1.5d0
        E(9, 3) = -dsqrt(2.d0)
        ! Matrix F
        F(1, 2) = 1.d0 / dsqrt(2.d0)
        F(2, 3) = 1.d0 / dsqrt(12.d0)
        F(3, 4) = -1.d0 / dsqrt(12.d0)
        F(4, 5) = -1.d0 / dsqrt(2.d0)
        F(6, 7) = -0.5d0
        F(7, 8) = 0.5d0
        F(1, 6) = -1.d0 / dsqrt(2.d0)
        F(2, 7) = 0.5d0
        F(3, 8) = dsqrt(3.d0) / 2.d0
        F(2, 9) = -dsqrt(2.d0 / 3.d0)
        F(6, 3) = dsqrt(3.d0) / 2.d0
        F(7, 4) = 0.5d0
        F(8, 5) = -1.d0 / dsqrt(2.d0)
        F(9, 4) = dsqrt(2.d0 / 3.d0)
        ! Matrix G
        G(1, 3) = -1.d0 / dsqrt(6.d0)
        G(2, 4) = -0.5d0
        G(3, 5) = -1.d0 / dsqrt(6.d0)
        G(1, 7) = 1.d0 / dsqrt(2.d0)
        G(2, 8) = 0.5d0
        G(1, 9) = -1.d0 / dsqrt(3.d0)
        G(6, 4) = -0.5d0
        G(7, 5) = -1.d0 / dsqrt(2.d0)
        G(6, 8) = 0.5d0
        G(9, 5) = -1.d0 / dsqrt(3.d0)
        ! Matrix Q20
        do i = 1, 9
            do j = 1, 9
                Q20(i, j) = Q0O * E(i, j) / 4.d0
            enddo
        enddo
        ! Matrix Q21
        do i = 1, 9
            do j = 1, 9
                Q21(i, j) = 3.d0 * Q0O * F(i, j) / (2.d0 * dsqrt(2.d0))
            enddo
        enddo
        ! Matrix Q22
        do i = 1, 9
            do j = 1, 9
                Q22(i, j) = 3.d0 * Q0O * G(i, j)
            enddo
        enddo
        ! Matrix Q2m1
        do i = 1, 9
            do j = 1, 9
                Q2m1(i, j) = (-1.d0) * Q21(j, i)
            enddo
        enddo
        ! Matrix Q2m2
        do i = 1, 9
            do j = 1, 9
                Q2m2(i, j) = Q22(j, i)
            enddo
        enddo
        ! Matrix A
        do i = 1, 9
            do j = 1, 9
                a1 = q10O2 / (2.d0 * rr**4) * (6.d0 * cc * Q20(i, j) + s * (Q21(i, j)&
                        - Q2m1(i, j)))
                a2 = 2.d0 * q20O2 / (48.d0 * rr**5) * (72.d0 * (3.d0 * cc * cc - 1.d0) * Q20(i, j)&
                        + 48.d0 * s * cc * (Q21(i, j) - Q2m1(i, j)) + 3.d0 * s * s * (Q22(i, j)&
                        + Q2m2(i, j)))
                !      print '(2i4,3g16.8)',i,j,a1,a2,a1+a2
                A(i, j) = a1 + a2
            enddo
        enddo
        do i = 1, 9
            do j = 1, 9
                vm(i, j) = A(i, j)
            enddo
        enddo
        do i = 1, Nd
            do j = 1, Nd
                if(dabs(vm(i, j) - E1(i, j))>1.e-8)then
                    !      write(6,*)'erreur de calcul de E1'
                    !      write(6,*)i,j,vm(i,j),E1(i,j)
                endif
                vmfso(i, j) = EFS(i, j) + (vm(i, j) + E2(i, j)) * autocm
                vmf(i, j) = (vm(i, j) + E2(i, j)) * autocm
            enddo
        enddo
        return
    end

    subroutine init(roo, cnee1, c6jj1, coef, ll)
        use splin

        implicit none

        integer :: ic, ipp
        real(kind = 8) :: roo, rooc, cnee1(nfj, nfj, 0:lmax, 0:2 * lmax + 1), c6jj1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)

        double precision :: coef(splin_nlt, splin_nr), rro(splin_nr), qq2(splin_nr)
        integer :: ll(splin_nlt, 6)

        !        parameter(nfj = 9, lmax = 2)
        !        common/cne/cnee1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        !        common/c6j/c6jj1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        !        common/dcc/ic
        !        common/droo/rooc

        ic = 0
        rooc = roo
        ! initialisation des matrices Cn en couplage jj
        ipp = 1
        call cnjel(roo, ic, ipp, cnee1)
        ! initialisation des matrices C6 en couplage jj
        ipp = 1
        call c6jmat(roo, ic, ipp, c6jj1)
        return
    end

    subroutine vpot(roo, rr, gam, vx, vxso)
        use splin

        implicit none

        integer :: Nd, ivec

        real(kind = 8) :: roo, rr, gam, vx, vxso
        real(kind = 8) :: VM(Ld, Ld), VL(Ld, Ld), VMSO(Ld, Ld)
        real(kind = 8) :: cnee1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)
        real(kind = 8) :: c6jj1(nfj, nfj, 0:lmax, 0:2 * lmax + 1)

        ivec = 1

        call buildmatrix(roo, rr, gam, VM, VMSO, cnee1, c6jj1)
        call hdiag(VM, Nd, ivec, Vl)
        vx = vm(9, 9)
        !       write(16,'(21g16.8)')roo,rr,gam,(vm(i,i),i=1,Nd)
        call hdiag(VMSO, Nd, ivec, Vl)
        vxso = vmso(9, 9)
        !       write(17,'(21g16.8)')roo, rr,gam,(vmso(i,i),i=1,Nd)
        return
    end

    subroutine matd(il, imm, imb, gm, dd)

        implicit none

        integer :: il, imm, imb, ii, mb, mm

        real(kind = 8) :: dcs, dsn, gm, cc, dd

        dcs = dcos(gm)
        dsn = dsin(gm)
        mb = imb - il
        mm = imm - il
        cc = 1.d0
        if(mb>mm)then
            cc = (-1)**(mb - mm)
            ii = mm
            mm = mb
            mb = ii
        endif
        if((iabs(mb)>iabs(mm)).or.(mb<0.and.mm<0))then
            ii = mm
            mm = -mb
            mb = -ii
        endif
        if(il==0)then
            dd = 1.d0
        else if(il==1)then
            if(mm==0.and.mb==0)then
                dd = dcos(gm)
            else if(mm==1.and.mb==-1)then
                dd = (1.d0 - dcos(gm)) / 2.d0
            else if(mm==1.and.mb==0)then
                dd = (-1.d0) * dsin(gm) / dsqrt(2.d0)
            else if(mm==1.and.mb==1)then
                dd = (1.d0 + dcos(gm)) / 2.d0
            else
                print *, 'erreur de mm,mb', il, mm, mb
                stop
            endif
        else if(il==2)then
            if(mm==0.and.mb==0)then
                dd = 1.5d0 * dcos(gm) * dcos(gm) - 0.5d0
            else if(mm==1.and.(mb + 1)==0)then
                dd = (1.d0 - dcos(gm)) * (2.d0 * dcos(gm) + 1.d0) / 2.d0
            else if(mm==1.and.mb==0)then
                dd = (-1.d0) * dsqrt(1.5d0) * dsin(gm) * dcos(gm)
            else if(mm==1.and.mb==1)then
                dd = (1.d0 + dcos(gm)) * (2.d0 * dcos(gm) - 1.d0) / 2.d0
            else if(mm==2.and.(mb + 2)==0)then
                dd = (1.d0 - dcos(gm))**2 / 4.d0
            else if(mm==2.and.(mb + 1)==0)then
                dd = (-1.d0) * (1.d0 - dcos(gm)) * dsin(gm) / 2.d0
            else if(mm==2.and.mb==0)then
                dd = dsqrt(6.d0) * dsin(gm) * dsin(gm) / 4.d0
            else if(mm==2.and.mb==1)then
                dd = (-1.d0) * (1.d0 + dcos(gm)) * dsin(gm) / 2.d0
            else if(mm==2.and.mb==2)then
                dd = (1.d0 + dcos(gm))**2 / 4.d0
            else
                print *, 'erreur de mm,mb', il, mm, mb
                stop
            endif
        else
            print *, 'erreur de il'
            stop
        endif
        dd = dd * cc
        return
    end

    SUBROUTINE HDIAG(H, N, IVEC, U)
        IMPLICIT None

        integer :: I, NM, N, NT, IVEC, J, IPU, IX, IXP, JX, &
                IC, ICV, ICW, IP, IXV, JXV, K
        real(kind = 8) :: H(Ld, Ld), U(Ld, Ld), DIA(Ld)
        real(kind = 8) :: ABE, NBS, ABH, COSN, COSV, DH12, DIFR, &
                H1, H12, H2, HBS, HGA, HGAJ, SINUS, TANG, TH, UGA, UT
        !     CALL ERRSET(208,300,-1)
        !     DIAGONALISATION EN DOUBLE PRECISION
        NT = 0
        NM = N - 1
        U(N, N) = 1.0
        IF(NM <= 0) then
            goto 99
        else
            goto 12
        END IF
        12 IF(IVEC == 0) then
            goto 10
        else
            goto 161
        END IF
        10 DO I = 1, NM
            U(I, I) = 1.0
            IPU = I + 1
            DO J = IPU, N
                U(I, J) = 0.
                U(J, I) = 0.
            end do
        end do
        161  DO IX = 1, N
            DIA(IX) = H(IX, IX)
        end do
        DO IX = 1, NM
            IXP = IX + 1
            DO JX = IXP, N
                ABH = DABS(H(IX, JX))
                104 H12 = H(IX, JX)
                H(IX, JX) = 0.
                H1 = H(IX, IX)
                H2 = H(JX, JX)
                DIFR = H1 - H2
                IF(IVEC == 0) then
                    go to 107
                else
                    go to 108
                END IF
                107 H(IX, IX) = H1 + H12
                H(JX, JX) = H1 - H12
                COSN = 0.707106781186547
                SINUS = COSN
                GO TO 109
                108 DH12 = H12 + H12
                TANG = DABS(DIFR) + DSQRT(DIFR * DIFR + DH12**2)
                TANG = DH12 / DSIGN(TANG, DIFR)
                COSV = 1. + TANG**2
                H(IX, IX) = (H1 + (DH12 + TANG * H2) * TANG) / COSV
                H(JX, JX) = (H2 - (DH12 - TANG * H1) * TANG) / COSV
                COSN = 1. / DSQRT(COSV)
                SINUS = TANG * COSN
                109 CONTINUE
                DO IC = 1, N
                    IXV = IX
                    ICV = IC
                    JXV = JX
                    ICW = IC
                    IF(IC - IX < 0) then
                        go to 204
                    elseif (IC - IX == 0) then
                        go to 103
                    else
                        go to 205
                    END IF
                    204 IXV = IC
                    ICV = IX
                    205 HGA = H(IXV, ICV)
                    IF(IC - JX < 0) then
                        go to 202
                    elseif (IC - JX == 0) then
                        go to 103
                    else
                        go to 106
                    END IF
                    202 ICW = JX
                    JXV = IC
                    106 HGAJ = H(JXV, ICW)
                    H(IXV, ICV) = HGA * COSN + HGAJ * SINUS
                    H(JXV, ICW) = HGAJ * COSN - HGA * SINUS
                    103 CONTINUE
                    IF(IVEC == 0) then
                        go to 210
                    else
                        go to 1021
                    END IF
                    210 UGA = U(IC, IX)
                    U(IC, IX) = SINUS * U(IC, JX) + COSN * UGA
                    U(IC, JX) = -SINUS * UGA + COSN * U(IC, JX)
                    1021 CONTINUE
                end do
            end do
        end do
        DO IX = 1, N
            IF(DABS(H(IX, IX) - DIA(IX)) - 1.E-10 <= 0) then
                go to 13
            else
                go to 170
            END IF
            13 CONTINUE
        end do
        GOTO 191
        170 NT = NT + 1
        IF(NT - 25 <= 0) then
            go to 161
        else
            go to 19
        END IF
        19 CONTINUE
        191  DO I = 1, NM
            IP = I + 1
            DO J = IP, N
                IF(H(I, I) - H(J, J) < 0) then
                    go to 212
                else
                    go to 209
                END IF
                212 TH = H(I, I)
                H(I, I) = H(J, J)
                H(J, J) = TH
                IF(IVEC == 0) then
                    go to 218
                else
                    go to 209
                END IF
                218 DO K = 1, N
                    UT = U(K, J)
                    U(K, J) = U(K, I)
                    U(K, I) = UT
                end do
                209 CONTINUE
            end do
        end do
        DO IX = 1, N
            HBS = DABS(H(IX, IX))
            NBS = HBS + 6.E-10
            ABE = NBS
            IF(DABS(ABE - HBS) - 6.E-10 <= 0) then
                go to 15
            else
                go to 16
            END IF
            15 H(IX, IX) = DSIGN(ABE, H(IX, IX))
            16 DO JX = 1, N
                HBS = DABS(U(IX, JX))
                NBS = HBS + 4.D-10
                ABE = NBS
                IF(DABS(ABE - HBS) - 5.E-10 <= 0) then
                    go to 18
                else
                    go to 17
                END IF
                18 U(IX, JX) = DSIGN(ABE, U(IX, JX))
                17 CONTINUE
            end do
        end do
        99 RETURN
    END

end module poteval