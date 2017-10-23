!Write by Liao Yuanda
!
!###############################################################

module modu_skx
contains
    subroutine calcu_Hf( L, N, M, Hf, t)
        !###############################################
        !(i,j) to (i, j_new); (i,j) to (i_new,j); N = i+(j-1)*L
        !###############################################
        implicit none
        integer :: L, N, M, i, j
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                p = i*L+j+1
                call triangular_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )

                Hf(p, p_i_pl) = -t
                Hf(p, p_i_mi) = -t
                Hf(p, p_j_pl) = -t
                Hf(p, p_i_mi) = -t
                Hf(p, p_i_pl_j_mi) = -t
                Hf(p, p_i_mi_j_pl) = -t

                Hf(p+N, p_i_pl+N) = -t
                Hf(p+N, p_i_mi+N) = -t
                Hf(p+N, p_j_pl+N) = -t
                Hf(p+N, p_i_mi+N) = -t
                Hf(p+N, p_i_pl_j_mi+N) = -t
                Hf(p+N, p_i_mi_j_pl+N) = -t
            end do
        end do
    end subroutine calcu_Hf

    subroutine calcu_H( N, M, Hf, H, J_k, fai, theta )
        implicit none
        integer :: N, M, i, j
        real*8 :: J_k, fai(N), theta(N), Hf(M,M)
        complex*16 :: H(M, M)
        do i = 1, M
            do j = 1, M
                H(i,j) = (0.0, 0.0)
            end do
        end do

        do i = 1, N 
            H(i, i) = -J_k*cos(theta(i)) 
            H(i+N, i+N) = J_k*cos(theta(i)) 
            H(i+N, i) = -J_k*sin(theta(i))*exp( cmplx(0.0, fai(i)) ) 
            H(i, i+N) = -J_k*sin(theta(i))*exp( cmplx(0.0, fai(i)*(-1)) ) 
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+H(i,j)
            end do
        end do

    end subroutine calcu_H

    subroutine calcu_Hp( p, L, J_H, lambda, fai, theta, Hp)
        implicit none
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L
        real*8 :: J_H, lambda, fai(L*L), theta(L*L), Hp
        Hp = 0.0
        call triangular_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )

        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_i_pl))*cos(fai(p)-fai(p_i_pl))+cos(theta(p))*cos(theta(p_i_pl)) )
        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_i_mi))*cos(fai(p)-fai(p_i_mi))+cos(theta(p))*cos(theta(p_i_mi)) )
        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_j_pl))*cos(fai(p)-fai(p_j_pl))+cos(theta(p))*cos(theta(p_j_pl)) )
        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_j_mi))*cos(fai(p)-fai(p_j_mi))+cos(theta(p))*cos(theta(p_j_mi)) )
        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_i_pl_j_mi))*cos(fai(p)-fai(p_i_pl_j_mi))+cos(theta(p))*cos(theta(p_i_pl_j_mi)) )
        Hp = Hp-J_H*( sin(theta(p))*sin(theta(p_i_mi_j_pl))*cos(fai(p)-fai(p_i_mi_j_pl))+cos(theta(p))*cos(theta(p_i_mi_j_pl)) )

        Hp = Hp-lambda*( sin(theta(p))*sin(fai(p))*cos(theta(p_i_pl))-cos(theta(p))*sin(theta(p_i_pl))*sin(fai(p_i_pl)) )
        Hp = Hp+lambda*( sin(theta(p))*sin(fai(p))*cos(theta(p_i_mi))-cos(theta(p))*sin(theta(p_i_mi))*sin(fai(p_i_mi)) )
        Hp = Hp-(sqrt(3.0)/2)*lambda*( cos(theta(p))*sin(theta(p_j_pl))*cos(fai(p_j_pl))&
            -sin(theta(p))*cos(fai(p))*cos(theta(p_j_pl)) )&
            -(1.0/2)*lambda*( sin(theta(p))*sin(fai(p))*cos(theta(p_j_pl))-cos(theta(p))*sin(theta(p_j_pl))*sin(fai(p_j_pl)) )
        Hp = Hp+(sqrt(3.0)/2)*lambda*( cos(theta(p))*sin(theta(p_j_mi))*cos(fai(p_j_mi))&
            -sin(theta(p))*cos(fai(p))*cos(theta(p_j_mi)) )&
            +(1.0/2)*lambda*( sin(theta(p))*sin(fai(p))*cos(theta(p_j_mi))-cos(theta(p))*sin(theta(p_j_mi))*sin(fai(p_j_mi)) )
        Hp = Hp-lambda*(1.0/2)*( sin(theta(p))*sin(fai(p))*cos(theta(p_i_pl_j_mi))&
            -cos(theta(p))*sin(theta(p_i_pl_j_mi))*sin(fai(p_i_pl_j_mi)) )&
            +lambda*(sqrt(3.0)/2)*( cos(theta(p))*sin(theta(p_i_pl_j_mi))*cos(fai(p_i_pl_j_mi))&
            -sin(theta(p))*cos(fai(p))*cos(theta(p_i_pl_j_mi)) )
        Hp = Hp+lambda*(1.0/2)*( sin(theta(p))*sin(fai(p))*cos(theta(p_i_mi_j_pl))&
            -cos(theta(p))*sin(theta(p_i_mi_j_pl))*sin(fai(p_i_mi_j_pl)) )&
            -lambda*(sqrt(3.0)/2)*( cos(theta(p))*sin(theta(p_i_mi_j_pl))*cos(fai(p_i_mi_j_pl))&
            -sin(theta(p))*cos(fai(p))*cos(theta(p_i_mi_j_pl)) )
        !write(*,*) Hp

    end subroutine calcu_Hp

    subroutine calcu_V( L, J_H, lambda, fai, theta, H, V )
        implicit none
        integer :: p, L, i, j
        real*8 :: J_H, lambda, fai(L*L), theta(L*L), H(L*L), V
        V = 0.0

        do i = 1,L
            do j = 1,L
                p = (i-1)*L+j
                call calcu_Hp( p, L, J_H, lambda, fai, theta, H(p))
                V = V+H(p)
            end do
        end do!
        V = V/2
        !write(*,*) V

    end subroutine calcu_V

    subroutine calcu_deltaV(L, N, B_z, theta, deltaV)
        implicit none
        integer :: L, i, j, N_ij, N
        real*8 :: B_z, theta(N), deltaV
        deltaV = 0
        do i = 0, L-1, 1
            do j = 0, L-1, 1
                N_ij = i*L+j+1
                deltaV = deltaV+B_z*cos(theta(N_ij))
            end do
        end do

    end subroutine calcu_deltaV

    subroutine calcu_W(M, H, beta, ln_Weight, Eig, Work, Rwork, miu )
        implicit none
        integer :: M, info, i, j
        real*8 :: beta, ln_Weight, miu
        complex*16 :: H(M,M), Work(3*M)
		complex*16 :: A(M, M)
        real*8 :: Eig(M),  Rwork(3*M-2)
        call zheev('N', 'U', M, H, M, Eig, Work, 3*M, Rwork, info)
        if (info /= 0) then
            write(*,*) 'info_erro =', info
        end if
        ln_Weight = 0.0 
        !write(*,*) Eig(:)
        do i = 1, M
            ln_Weight = ln_Weight+log( 1+exp( (-1)*beta*(Eig(i)-miu) ) )
        end do
    end subroutine calcu_W


end module modu_skx

!function to less dimension
integer function di(d)
    integer :: d
    if (d == 1) then
        di = 2
    end if
    if ( d == 2) then
        di =1
    end if
end function


program skyrmionlattice
    use modu_skx 
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_k, the parameter if H
    !##########################################
    integer, parameter :: L = 20, N = L*L, M = 2*N, num_T_para = 20, num_MC = 4*10**4
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, r, d, di, dd, ij, aver_num
    integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl
    real*8 :: ran, J_k, beta, miu, Tc, x, Weight_fermi, Weight_fermi_new, a, b, t1, t2
    real*8 :: V, deltaV, ln_Weight_fermi_new, J_H, lambda, B_z, pa, Hp, energy_aver, energy_aver_B
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(num_T_para) :: T_para
    real*8, dimension(N) :: Hc
    real*8, dimension(M) :: Eig
    real*8, dimension(3*M-2) :: Rwork
    real*8, dimension(M,M) :: Hf
    complex*16, dimension(3*M) :: Work
    complex*16, dimension(M,M) :: H 
    character(len=2) :: cha
    open( unit =12, file = 'miu.in')
    read(12,*) miu
    J_k = t
    J_H = t
    lambda = 0.5*PI*t
    !open(unit = 12, file = 'B_z.in')
    !read(12,*) pa
    !close(12)
    !B_z = pa*J_H
    B_z = 0.8*J_H

    d = 1
    aver_num = 100
    T_para = (/ 1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05/)

    open(unit = 15, file = 'T_energyaver.out')
	
	!#############
	call cpu_time(t1)
	write(*,*) 'Begin at :', t1
	!##################
 
   call initran(1)

    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI
    end do

    call calcu_V( L, J_H, lambda, fai(d,:), theta(d,:), Hc, V )
 
    call initran(1)
    call calcu_Hf( L, N, M, Hf, t )

do ii = 1, num_T_para
    energy_aver = 0.0
    energy_aver_B = 0.0
    write(cha,'(i2)') ii
    open(unit = 13, file = 'energy'//adjustl(trim(cha))//'.out' )
    open( unit = 14, file = 'lattice'//adjustl(trim(cha))//'.out' )
    write(*,*) num_T_para-ii
    Tc = T_para(ii)*J_H
    beta = 1/Tc
	call calcu_H( N, M, Hf, H, J_k, fai(d,:), theta(d,:))
    call calcu_W(M, H, beta, Weight_fermi, Eig, Work, Rwork, miu )

    do i = 1, num_MC

        dd = di(d)
        x = ran()
        r = int(x*N)+1
        fai(dd,:) = fai(d,:)
        theta(dd,:) = theta(d,:)
        x = ran()
        fai(dd, r) = fai(d,r)+(x-0.5)*2*2*PI/12
        if ( fai(dd,r) < 0 ) then
            fai(dd,r) = fai(dd,r)+2*PI
        end if
        if ( fai(dd,r) > 2*PI ) then
            fai(dd,r) = fai(dd,r)-2*PI
        end if
        x = ran()
        theta(dd, r) = theta(d,r)+(x-0.5)*2*PI/6
        if ( theta(dd,r) < 0 ) then
            theta(dd,r) = theta(dd,r)+PI
        end if
        if (theta(dd,r) > PI ) then
            theta(dd,r) = theta(dd,r)-PI
        end if
		call calcu_H( N, M, Hf, H, J_k, fai(dd,:), theta(dd,:))
        call calcu_W(M, H, beta, Weight_fermi_new, Eig, Work, Rwork, miu )
        call calcu_Hp( r, L, J_H, lambda, fai(dd,:), theta(dd,:), Hp)
        deltaV = B_z*( cos(theta(dd,r))-cos(theta(d,r)) )
        x = ran()

        !write(*,*) x, '<?', exp( Weight_fermi_new-Weight_fermi-beta*(Hp-Hc(r)+deltav) )
        !write(*,*) Weight_fermi, Weight_fermi_new
        !write(*,*) V, V+Hp-Hc(r)
        if ( x > exp( Weight_fermi_new-Weight_fermi-beta*(Hp-Hc(r)+deltaV) ) ) then
            fai(dd,r) = fai(d,r)
            theta(dd, r) = theta(d, r)
        else
        	Weight_fermi = Weight_fermi_new
            V = V+Hp-Hc(r)
            Hc(r) = Hp
            call triangular_nearest( r, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )
            call calcu_Hp( p_i_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_i_pl))
            call calcu_Hp( p_i_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_i_mi))
            call calcu_Hp( p_j_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_j_pl))
            call calcu_Hp( p_j_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_j_mi))
            call calcu_Hp( p_i_pl_j_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_i_pl_j_mi))
            call calcu_Hp( p_i_mi_j_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), Hc(p_i_mi_j_pl))
        end if

        d = dd
        if ( i > num_MC-aver_num*N ) then
            if ( mod(i+aver_num*N-num_MC, N) == 0 ) then
                energy_aver = energy_aver+V
                call calcu_deltaV(L, N, B_z, theta(d,:), deltaV)
                energy_aver_B = energy_aver_B+V+deltav
                do j = 1, N
                    write(14,*) theta(d,j), fai(d,j)
                end do
            end if
        end if
        write(13,*) V
    end do

    close(13)
    close(14)

    energy_aver = energy_aver/aver_num
    energy_aver_B = energy_aver_B/aver_num
    write(15,*) T_para(ii), energy_aver, energy_aver_B
end do

    close(15)
	!#############
	call cpu_time(t2)
	write(*,*) 'End at :', t2
	!##################
    write(*,*) 'Timecost:', t2-t1
end program skyrmionlattice
