!Write by Liao Yuanda
!
!###############################################################

module dou_ex
contains

    subroutine calcu_V( L, N, J_H, lambda, fai, theta, V )
        implicit none
        integer :: L, N, i, j, i_minus, i_plus, j_minus, j_plus, N_i_minus, N_i_plus, N_j_minus, N_j_plus, N_ij, N_i_plus_j_minus, N_i_minus_j_plus
        real*8 :: J_H, lambda, fai(N), theta(N), V
        V = 0
        do i = 0, L-1, 1
            do j = 0, L-1, 1
                !period of boundary 
                i_plus = mod(i+1, L)
                j_plus = mod(j+1, L)
                if (i==0) then
                    i_minus = L-1
                else
                    i_minus = i-1
                end if 
                if (j==0) then
                    j_minus = L-1
                else
                    j_minus = j-1
                end if 
                N_ij = i*L+j+1
                N_i_plus = i_plus*L+j+1
                N_j_plus = i*L+j_plus+1
                N_i_plus_j_minus = i_plus*L+j_minus+1

                V = V-J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_plus))*cos(fai(N_i_plus))
                V = V-J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_j_plus))*cos(fai(N_j_plus))
                V = V-J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_plus_j_minus))*cos(fai(N_i_plus_j_minus))

                V = V-J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_plus))*sin(fai(N_i_plus))
                V = V-J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_j_plus))*sin(fai(N_j_plus))
                V = V-J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus))

                V = V-J_H*cos(theta(N_ij))*cos(theta(N_i_plus))
                V = V-J_H*cos(theta(N_ij))*cos(theta(N_j_plus))
                V = V-J_H*cos(theta(N_ij))*cos(theta(N_i_plus_j_minus))

                V = V+lambda*( sin(theta(N_ij))*sin(fai(N_ij))*cos(theta(N_i_plus))-cos(theta(N_ij))&
                    *sin(theta(N_i_plus))*sin(fai(N_i_plus)) )
                V = V+lambda*( (1.0/2)*( sin(theta(N_ij))*sin(fai(N_ij))*cos(theta(N_j_plus))&
                    -cos(theta(N_ij))*sin(theta(N_j_plus))*sin(fai(N_j_plus)) ) &
                    + sqrt(3.0)/2*( sin(theta(N_ij))*cos(fai(N_ij))*cos(theta(N_j_plus))-cos(theta(N_ij))*sin(theta(N_j_plus))*cos(fai(N_j_plus)) ) )
                V = V+lambda*( sin(theta(N_ij))*sin(fai(N_ij))*cos(theta(N_i_plus_j_minus))-cos(theta(N_ij))&
                    *sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus)) )&
                    +lambda*(-1.0)*( (1.0/2)*( sin(theta(N_ij))*sin(fai(N_ij))*cos(theta(N_i_plus_j_minus))&
                    -cos(theta(N_ij))*sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus)) ) &
                    + sqrt(3.0)/2*( sin(theta(N_ij))*cos(fai(N_ij))*cos(theta(N_i_plus_j_minus))&
                    -cos(theta(N_ij))*sin(theta(N_i_plus_j_minus))*cos(fai(N_i_plus_j_minus)) ) )

            end do
        end do


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
    
end module dou_ex

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
    use dou_ex 
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_k, the parameter if H
    !##########################################
    integer, parameter :: L = 40, N = L*L, num_MC = 5*10**5, num_T_para = 21
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, p, d, di, dd, ij, aver_num
    real*8 :: ran, J_k, beta, miu, Tc, x, energy
    real*8 :: V, deltaV, J_H, lambda, B_z, pa, V_new, deltaV_new
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(num_T_para) :: T_para
    character(len=2) :: cha
    !open( unit = 12, file = 'para.in' )
    !read(12,*) pa
    !pa = 55.0/89
    d = 1
    lambda = 0.4*PI
    J_H = lambda/11.6239
    !B_z = 0.001*t
    !deltaV = 0
    !V = 0
    aver_num = 100

    T_para = (/1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2 ,0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01,&
        0.006, 0.002/)
    call initran(1)
    open( unit = 13, file = 'lattice.out' )

    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI
    end do
    !ij = 0
do ii = 1, num_T_para 
    write(cha,'(i2)') ii
    open(unit = ii+20, file = 'energy'//adjustl(trim(cha))//'.out' ) 
    write(*,*) num_T_para-ii
    Tc = T_para(ii)*J_H
    beta = 1/Tc
    p = 0
    call calcu_V( L, N, J_H, lambda, fai(d,:), theta(d,:), V )
    !call calcu_deltaV(L, N, B_z, theta(d,:), deltaV)
    !ln_Weight = -beta*V

    do i = 1, num_MC 
        dd = di(d)
        !x = ran()
        !p = int(x*N)+1 
        p = p+1
        !Write(*,*) p
        fai(dd,:) = fai(d,:)
        theta(dd,:) = theta(d,:)
        x = ran()
        fai(dd, p) = fai(d,p)+(x-0.5)*2*2*PI/12
        if ( fai(dd,p) < 0 ) then
            fai(dd,p) = fai(dd,p)+2*PI
        end if
        if ( fai(dd,p) > 2*PI ) then
            fai(dd,p) = fai(dd,p)-2*PI
        end if

        x = ran()
        theta(dd, p) = theta(d,p)+(x-0.5)*2*PI/6
        if ( theta(dd,p) < 0 ) then
            theta(dd,p) = theta(dd,p)+PI
        end if
        if (theta(dd,p) > PI ) then
            theta(dd,p) = theta(dd,p)-PI
        end if
        call calcu_V( L, N, J_H, lambda, fai(dd,:), theta(dd,:), V_new )
        !call calcu_deltaV(L, N, B_z, theta(dd,:), deltaV)
        !ln_Weight_new = -beta*(V+deltaV)

        x = ran()
        !write(*,*) x, '<?', exp( -beta*V_new+beta*V )
        !write(*,*) -beta*V, -beta*V_new
        if ( x > exp( -beta*V_new+beta*V ) ) then
            fai(dd,p) = fai(d,p)
            theta(dd, p) = theta(d, p)
        else
            !ln_Weight = ln_Weight_new
            V = V_new
        end if
        d = dd
        !write(*,*) i/(num_MC/100)
        if (ii == num_T_para .AND. i > num_MC-aver_num*N ) then
            if ( mod(i+aver_num*N-num_MC, N) == 0 ) then
                !energy = energy+V
                !ij = ij+1
                !write(13,*) ij
                do j = 1, N
                    write(13,*) j, theta(d,j), fai(d,j)
                end do
            end if
        end if
        if ( p == N ) then
            p = 0
        end if
        write(ii+20,*) V

    end do
    close(ii+20)
end do
!energy = energy/aver_num
!open(18,file = 'energy.out')
!write(18,*) energy


    write(*,*) 'end'
end program skyrmionlattice
