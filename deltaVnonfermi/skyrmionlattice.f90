!Write by Liao Yuanda
!
!###############################################################

module modu
contains
    subroutine calcu_p_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )
        implicit none
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L
        integer :: i, j, i_pl, j_pl, i_mi, j_mi

        i = int((p-1)/L)
        j = mod(p-1,L)
        i_pl = mod(i+1,L)
        j_pl = mod(j+1,L)
        if (i==0) then
            i_mi = L-1
        else
            i_mi = i-1
        end if
        if (j==0) then
            j_mi = L-1
        else
            j_mi = j-1
        end if

        p_i_pl = i_pl*L+j+1
        p_i_mi = i_mi*L+j+1
        p_j_pl = i*L+j_pl+1
        p_j_mi = i*L+j_mi+1
        p_i_pl_j_mi = i_pl*L+j_mi+1
        p_i_mi_j_pl = i_mi*L+j_pl+1

    end subroutine calcu_p_nearest 

    subroutine calcu_Hp( p, L, J_H, lambda, fai, theta, Hp)
        implicit none
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L
        real*8 :: J_H, lambda, fai(L*L), theta(L*L), Hp
        Hp = 0.0
        call calcu_p_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )

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
    
end module modu

real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
    implicit none

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64

    ran64=ran64*mul64+add64
    ran=0.5d0+dmu64*dble(ran64)

end function ran
!----------------!

!---------------------!
subroutine initran(w)
!---------------------!
    implicit none

    integer(8) :: irmax
    integer(4) :: w,nb,b

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64
                        
    irmax=2_8**31
    irmax=2*(irmax**2-1)+1
    mul64=2862933555777941757_8
    add64=1013904243
    dmu64=0.5d0/dble(irmax)

    open(10,file='seed.in',status='old')

    read(10,*)ran64
    close(10)
    if (w.ne.0) then
        open(10,file='seed.in',status='unknown')
        write(10,*) abs((ran64*mul64)/5+5265361)
        close(10)
    end if

end subroutine initran
!----------------------!

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
    use  modu
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_k, the parameter if H
    !##########################################
    integer, parameter :: L = 40, N = L*L, num_MC = 4*10**5, num_T_para = 30
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, r, d, di, dd, ij, aver_num
    integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl
    real*8 :: ran, beta, Tc, x, Hp, B_z
    real*8 :: V, deltaV, J_H, lambda, pa, energy_aver, energy_aver_B
    integer, dimension(1) :: mloc
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(num_T_para) :: T_para
    real*8, dimension(N) :: H
    character(len=2) :: cha
    d = 1
    lambda = 0.4*PI*t
    J_H = t
    aver_num = 100
    open(unit = 12, file = 'B_z.in')
    read(12,*) pa
    close(12)
    B_z = pa*t
    T_para = (/ 1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08,&
        0.06, 0.04, 0.02, 0.008, 0.006, 0.004, 0.002, 0.0008, 0.0006, 0.0004/)
    open(unit = 15, file = 'T_energyaver.out')
    call initran(1)

    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI
    end do
    call calcu_V( L, J_H, lambda, fai(d,:), theta(d,:), H, V )
do ii = 1, num_T_para
    energy_aver = 0.0
    energy_aver_B = 0.0
    write(cha,'(i2)') ii
    open(unit = 13, file = 'energy'//adjustl(trim(cha))//'.out' )
    open( unit = 14, file = 'lattice'//adjustl(trim(cha))//'.out' )
    write(*,*) num_T_para-ii
    Tc = T_para(ii)*J_H
    beta = 1/Tc

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

        call calcu_Hp( r, L, J_H, lambda, fai(dd,:), theta(dd,:), Hp)
        deltaV = B_z*( cos(theta(dd,r))-cos(theta(d,r)) )

        x = ran()
        !write(*,*) x, '<?', exp( -beta*(Hp-H(r)+deltav) )
        if ( x > exp( -beta*(Hp-H(r)+deltaV) ) ) then
            fai(dd,r) = fai(d,r)
            theta(dd, r) = theta(d, r)
        else
            V = V+Hp-H(r)
            H(r) = Hp
            call calcu_p_nearest( r, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )
            call calcu_Hp( p_i_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_i_pl))
            call calcu_Hp( p_i_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_i_mi))
            call calcu_Hp( p_j_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_j_pl))
            call calcu_Hp( p_j_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_j_mi))
            call calcu_Hp( p_i_pl_j_mi, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_i_pl_j_mi))
            call calcu_Hp( p_i_mi_j_pl, L, J_H, lambda, fai(dd,:), theta(dd,:), H(p_i_mi_j_pl))
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

    write(*,*) 'end'
end program skyrmionlattice
